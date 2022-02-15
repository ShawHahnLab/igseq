from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from igseq import record
from .util import TestBase

class TestRecordHandler(TestBase):

    def setUp(self):
        super().setUp()
        self.handler = record.RecordHandler(self.path/"example.fasta")

    def tearDown(self):
        self.handler.close()
        super().tearDown()

    def test_infer_fmt(self):
        self.assertEqual(self.handler.fmt, "fa")
        fmt = self.handler.infer_fmt("fq")

    def test_encode_record(self):
        obj = self.handler.encode_record({"sequence_id": "id", "sequence": "ACTG"})
        self.assertEqual(obj.id, "id")
        self.assertEqual(str(obj.seq), "ACTG")

    def test_decode_record(self):
        rec = self.handler.decode_record({"sequence_id": "id", "sequence": "ACTG"})
        self.assertEqual(rec["sequence_id"], "id")
        self.assertEqual(rec["sequence"], "ACTG")
        rec = self.handler.decode_record(SeqRecord(Seq("ACTG"), id="id"))
        self.assertEqual(rec["sequence_id"], "id")
        self.assertEqual(rec["sequence"], "ACTG")

    def test_encode_phred(self):
        txt = self.handler.encode_phred([33, 33, 33])
        self.assertEqual(txt, "BBB")

    def test_decode_phred(self):
        nums = self.handler.decode_phred("BBB")
        self.assertEqual(nums, [33, 33, 33])

    def test_open(self):
        # open() is not implemented
        with self.assertRaises(NotImplementedError):
            self.handler.open()

    def test_close(self):
        # does nothing, by default, as no handle is open
        self.handler.close()


class TestRecordReader(TestRecordHandler):

    def setUp(self):
        super().setUp()
        self.handler = record.RecordReader(self.path/"example.fasta")

    def test_open(self):
        self.handler.open()

    def test_close(self):
        # no effect if not yet open
        self.handler.close()
        # will close if open
        self.handler.open()
        self.handler.close()

    def test_reading(self):
        # always gives one dictionary for each record, whatever the input type
        self.handler.open()
        recs = []
        for rec in self.handler:
            recs.append(rec)
        self.assertEqual(recs, [{"sequence_id": "id", "sequence": "ACTG"}])
        self.handler.close()
        # needs to be opened first
        recs = []
        for rec in self.handler:
            recs.append(rec)
            self.assertEqual(recs, [])

    def test_context_manager(self):
        # As a context manager it will open and close automatically
        with self.handler:
            self.assertFalse(self.handler.handle.closed)
        self.assertTrue(self.handler.handle.closed)


class TestRecordWriter(TestRecordHandler):

    def setUp(self):
        super().setUp()
        self.handler = record.RecordWriter(self.tmp/"example.fasta")

    def test_open(self):
        self.handler.open()
        self.assertTrue((self.tmp/"example.fasta").exists())

    def test_write(self):
        self.handler.open()
        self.handler.write({"sequence_id": "id", "sequence": "ACTG"})
        self.handler.close()
        self.assertTxtsMatch(self.path/"example.fasta", self.tmp/"example.fasta")

    def test_context_manager(self):
        with self.handler:
            self.assertFalse(self.handler.handle.closed)
        self.assertTrue(self.handler.handle.closed)
        with self.handler:
            self.handler.write({"sequence_id": "id", "sequence": "ACTG"})
        self.assertTxtsMatch(self.path/"example.fasta", self.tmp/"example.fasta")

class TestRecordWriterHandle(TestRecordWriter):
    """Test RecordWriter with an existing file handle"""

    def setUp(self):
        super().setUp()
        self.fobj = open(self.tmp/"example.fasta", "wt")
        self.handler = record.RecordWriter(self.fobj)

    def tearDown(self):
        if not self.fobj.closed:
            self.fobj.close()
        super().tearDown()

    def test_context_manager(self):
        # this time, it shouldn't automatically close the file since it's
        # already open.
        with self.handler:
            self.assertFalse(self.handler.handle.closed)
        self.assertFalse(self.handler.handle.closed)
        # The context manager is still helpful because it'll automatically call
        # close which will automatically call flush (even though it won't
        # close)
        with self.handler:
            self.handler.write({"sequence_id": "id", "sequence": "ACTG"})
        self.assertTxtsMatch(self.path/"example.fasta", self.tmp/"example.fasta")
