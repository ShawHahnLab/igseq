"""
Special case of one-off external commands for misc tests.

Each calls a bash script that should exit 0 if it worked, nonzero otherwise.
"""

import inspect
from subprocess import run
from .util import TestBase

class TestCmds(TestBase):
    """Misc tests contained in standalone bash scripts."""

    def test_broken_pipe(self):
        """Test that broken pipes are handled gracefully for an igblast call."""
        self.runcmd()

    def runcmd(self):
        """Run the test bash script associated with the calling function.

        This is either clever or idiotic.  Time will tell.
        """
        frame = inspect.currentframe()
        parent = inspect.getouterframes(frame)[1][0]
        func = inspect.getframeinfo(parent).function
        script = self.path / (func.removeprefix("test_") + ".sh")
        run(script, check=True)
