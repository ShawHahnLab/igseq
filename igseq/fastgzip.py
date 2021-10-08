"""
gzip reading and writing via the external gzip programs.
"""

from subprocess import Popen, PIPE
import types

GZIP = "gzip"
GUNZIP = "gunzip"

def fast_gzip_open(path, mode="rt"):
    """open()-like wrapper for fast_gzip/fast_gunzip"""
    if "w" in mode:
        return fast_gzip(path, mode)
    return fast_gunzip(path, mode)

def fast_gunzip(path, mode="rt"):
    """Open a gzip file for reading via the external gunzip executable."""
    text = "b" not in mode
    with open(path, "rb") as f_in:
        # Here we'll add some extra magic to the returned file handle so that,
        # when its close method is called, we handle some cleanup with the
        # underlying file handle and the running process.  (Just using the
        # context manager feature of Popen here doesn't seem to work.)
        # pylint: disable=consider-using-with
        proc = Popen([GUNZIP], stdin=f_in, stdout=PIPE, text=text)
        proc.stdout.gz_handle = f_in
        proc.stdout.proc = proc
        def method_close(self):
            # close gunzip's stdin (the gz stream)
            self.gz_handle.close()
            # wait for gunzip to exit
            self.proc.wait()
            # close gunzip's stdout
            self.close_real()
        proc.stdout.close_real = proc.stdout.close
        proc.stdout.close = types.MethodType(method_close, proc.stdout)
        return proc.stdout

def fast_gzip(path, mode="wt"):
    """Open a gzip file for writing via the external gzip executable."""
    text = "b" not in mode
    with open(path, "wb") as f_out:
        # similar setup to fast_gunzip, just organized in reverse.
        # pylint: disable=consider-using-with
        proc = Popen([GZIP], stdin=PIPE, stdout=f_out, text=text)
        proc.stdin.gz_handle = f_out
        proc.stdin.gz_proc = proc
        def method_close(self):
            # close gzip's stdin for real
            self.close_real()
            # wait for gzip to exit
            self.gz_proc.wait()
            # close gzip's stdout
            self.gz_handle.close()
        proc.stdin.close_real = proc.stdin.close
        proc.stdin.close = types.MethodType(method_close, proc.stdin)
        return proc.stdin
