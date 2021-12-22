import logging
from .record import RecordReader, RecordWriter

def convert(input_fp, output_fp, input_fmt=None, output_fmt=None, dummyqual=None):
    with RecordReader(input_fp, input_fmt) as reader, \
        RecordWriter(output_fp, output_fmt, dummyqual=dummyqual) as writer:
        for record in reader:
            writer.write(record)
