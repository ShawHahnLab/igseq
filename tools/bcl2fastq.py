#!/usr/bin/env python

"""
Just sanity-checking some things for fastq2bcl.py.
"""

import re
import struct
from pathlib import Path

def parse_cluster_byte(raw):
    byte = struct.unpack("B", raw)[0]
    base = ["A", "C", "G", "T"][byte & 0x3]
    qual = chr((byte >> 2) + 33)
    return (base, qual)

def load_filter_data(path):
    with open(path, "rb") as f_in:
        filter_data = {}
        filter_data["prefix"] = f_in.read(4) # always zero
        filter_data["version"] = struct.unpack("<I", f_in.read(4))
        filter_data["cluster_count"] = struct.unpack("<I", f_in.read(4))
        cluster_filter_pass = []
        while True:
            byte = f_in.read(1)
            if not byte:
                break
            cluster_filter_pass.append(struct.unpack("B", byte)[0])
        filter_data["cluster_filter_pass"] = cluster_filter_pass
    return filter_data

def bcl2fastq(rundir):
    rundir = Path(rundir)
    cycle_dirs = [path for path in (rundir/"Data/Intensities/BaseCalls/L001").glob("*") if path.is_dir()]
    cycle_dirs = sorted(cycle_dirs, key=lambda p: int(re.sub(r"C([0-9]+)\.1", r"\1", str(p.name))))

    filter_data = load_filter_data(rundir/"Data/Intensities/BaseCalls/L001/s_1_1101.filter")

    handles = [open(c/"s_1_1101.bcl", "rb") for c in cycle_dirs]
    for handle in handles:
        # skip over the cluster number
        handle.read(4)
    idx = 0
    while True:
        try:
            pairs = [parse_cluster_byte(hndl.read(1)) for hndl in handles]
        except struct.error:
            break
        filter_status = filter_data["cluster_filter_pass"][idx]
        print(f"@read{idx} pass={filter_status}")
        print("".join([p[0] for p in pairs]))
        print("+")
        print("".join([p[1] for p in pairs]))
        idx += 1

if __name__ == "__main__":
    bcl2fastq(sys.argv[1])
