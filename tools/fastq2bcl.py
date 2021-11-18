#!/usr/bin/env python

"""
Mock up a bcl2fastq-able run directory from just a I1/R1/R2 fastq.gz trio.

Assumptions:

 * layout is R1+I1+R2
 * only lane 1
 * only tile 1101

https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf
"""

import re
import sys
import gzip
import struct
from pathlib import Path
from Bio import SeqIO

def parse_seqdesc_fields(txt):
    match = re.match(
        r"([A-Za-z0-9_]+):([0-9]+):([-A-Za-z0-9]*):([0-9]+):([0-9]+):([0-9]+):([0-9]+) "
        r"([12]):([YN]):([0-9]+):([A-Z]+)",
        txt)
    if not match:
        raise ValueError("ID not recognized: %s" % txt)
    keys = [
        "instrument", "run_number", "flowcell_id", "lane", "tile", "x_pos", "y_pos", "read",
        "is_filtered", "control_number", "barcode_sequence"]
    fields = {key: match.group(idx+1) for idx, key in enumerate(keys)}
    return fields

def mock_run_id(seqdesc_txt):
    fields = parse_seqdesc_fields(seqdesc_txt)
    run_id = "YYMMDD_" + \
        fields['instrument'] + "_" + \
        fields['run_number'].zfill(4) + "_" + \
        fields['flowcell_id']
    return run_id

def write_run_info_xml(xml_out, cycles_i1, cycles_r1, cycles_r2, desc):
    xml_out = Path(xml_out)
    xml_out.parent.mkdir(exist_ok=True, parents=True)
    fields = parse_seqdesc_fields(desc)
    run_id = mock_run_id(desc)
    run_number = fields["run_number"]
    flowcell_id = fields["flowcell_id"]
    instrument = fields["instrument"]
    with open(xml_out, "wt") as f_out:
        f_out.write('<?xml version="1.0"?> \n'
        '<RunInfo xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" Version="2">\n'
        f'  <Run Id="{run_id}" Number="{run_number}">\n'
        f'    <Flowcell>{flowcell_id}</Flowcell>\n'
        f'    <Instrument>{instrument}</Instrument>\n'
        f'    <Date>YYMMDD</Date>\n'
        f'    <Reads>\n'
        f'      <Read NumCycles="{cycles_r1}" Number="1" IsIndexedRead="N" />\n'
        f'      <Read NumCycles="{cycles_i1}" Number="2" IsIndexedRead="Y" />\n'
        f'      <Read NumCycles="{cycles_r2}" Number="3" IsIndexedRead="N" />\n'
        f'    </Reads>\n'
        '    <FlowcellLayout LaneCount="1" SurfaceCount="1" SwathCount="1" TileCount="1" />\n'
        '  </Run>\n'
        '</RunInfo>\n')

def write_filter(outdir, cluster_count):
    path = Path(outdir)/"Data/Intensities/BaseCalls/L001/s_1_1101.filter"
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wb") as f_out:
        f_out.write(bytes([0, 0, 0, 0])) # prefix, always zero
        f_out.write(bytes([3, 0, 0, 0])) # version, just one I saw
        f_out.write(struct.pack("<I", cluster_count))
        # every cluster passes!
        f_out.write(bytes([1]*cluster_count))

def write_control(outdir, cluster_count):
    path = Path(outdir)/"Data/Intensities/BaseCalls/L001/s_1_1101.control"
    path.parent.mkdir(exist_ok=True, parents=True)
    with open(path, "wb") as f_out:
        f_out.write(bytes([0, 0, 0, 0])) # "Zero value (for backwards compatibility)"
        f_out.write(bytes([2, 0, 0, 0])) # "Format version number"
        f_out.write(struct.pack("<I", cluster_count)) # "Number of clusters"
        # two bytes for each cluster
        f_out.write(bytes([0, 0]*cluster_count))

def write_bcls_and_stats(outdir, trios, cycles):
    joined = []
    for trio in trios:
        seqs = [str(rec.seq) for rec in trio]
        quals = [rec.letter_annotations["phred_quality"] for rec in trio]
        joined.append((seqs[1] + seqs[0] + seqs[2], quals[1] + quals[0] + quals[2]))
    # write cluster counts first
    for cycle in range(cycles):
        cycledir = outdir/f"Data/Intensities/BaseCalls/L001/C{cycle+1}.1"
        cycledir.mkdir(exist_ok=True, parents=True)
        with open(cycledir/"s_1_1101.bcl", "wb") as f_out:
            f_out.write(struct.pack("<I", len(trios)))
    # write individual bases across all clusters for each cycle
    for cycle in range(cycles):
        cycledir = outdir/f"Data/Intensities/BaseCalls/L001/C{cycle+1}.1"
        for basecalls, qualscores in joined:
            bcl_byte = encode_cluster_byte(basecalls[cycle], qualscores[cycle])
            with open(cycledir/"s_1_1101.bcl", "ab") as f_out:
                f_out.write(bcl_byte)
        with open(cycledir/"s_1_1101.stats", "wb") as f_out:
            # can I get away with this?
            f_out.write(bytes([0]*108))

def encode_cluster_byte(base, qual):
    if base == "N":
        return bytes([0]) # no call
    qual = qual << 2
    base = ["A", "C", "G", "T"].index(base)
    return bytes([qual | base])

def write_locs(outdir, trios):
    path = Path(outdir)/"Data/Intensities/L001/s_1_1101.locs"
    path.parent.mkdir(exist_ok=True, parents=True)
    #  looks like length is 12 prefix bytes and then 8 bytes for every cluster,
    # as a float for X and a float for Y.
    # example first 12 bytes from a run:
    # 01 00 00 00 = ?
    # 00 00 80 3f = ?
    # 2c 0c 06 00 = number of clusters in the tile
    # those first 8 are the same across locs files from other runs, sequencers,
    # etc. soooo I guess we'll just keep them
    with open(path, "wb") as f_out:
        f_out.write(bytes([1, 0, 0, 0, 0, 0, 0x80, 0x3f]))
        f_out.write(struct.pack("<I", len(trios)))
        for trio in trios:
            fields = parse_seqdesc_fields(trio[0].description)
            f_out.write(encode_loc_bytes(fields["x_pos"], fields["y_pos"]))

def encode_loc_bytes(x_pos, y_pos):
    x_bytes = struct.pack("<f", (int(x_pos) - 1000)/10)
    y_bytes = struct.pack("<f", (int(y_pos) - 1000)/10)
    return x_bytes + y_bytes

def fastq2bcl(fqgz_i1, fqgz_r1, fqgz_r2, outdir=None):
    trios = []
    with gzip.open(fqgz_i1, "rt") as f_i1, \
        gzip.open(fqgz_r1, "rt") as f_r1, \
        gzip.open(fqgz_r2, "rt") as f_r2:
        for trio in zip(SeqIO.parse(f_i1, "fastq"), \
            SeqIO.parse(f_r1, "fastq"), \
            SeqIO.parse(f_r2, "fastq")):
            if trio[0].id != trio[1].id or trio[1].id != trio[2].id:
                raise ValueError(
                    "Seq ID mismatch for I1/R1/R2: %s/%s/%s" % (trio[0].id, trio[1].id, trio[2].id))
            trios.append(trio)
    if not outdir:
        outdir = mock_run_id(trios[0][0].description)
    outdir = Path(outdir)
    cycles_i1 = len(trios[0][0])
    cycles_r1 = len(trios[0][1])
    cycles_r2 = len(trios[0][2])
    cycles = cycles_r1 + cycles_i1 + cycles_r2
    desc = trios[0][0].description
    write_run_info_xml(outdir/"RunInfo.xml", cycles_i1, cycles_r1, cycles_r2, desc)
    write_filter(outdir, cluster_count=len(trios))
    write_control(outdir, cluster_count=len(trios))
    write_locs(outdir, trios)
    write_bcls_and_stats(outdir, trios, cycles)

if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    fastq2bcl(*sys.argv[1:])
