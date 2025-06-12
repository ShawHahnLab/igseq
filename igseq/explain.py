"""
Show explanations for concepts used in multiple commands.

The explain command takes a key word (currently just "samples") and shows a
summary explanation.  Run without arguments to see a list of available entries.
"""

from .util import DATA, IgSeqError
from .show import show_text

def __load_explanations(path):
    paths = sorted((DATA/"explain").glob("*.txt"))
    paths = {path.stem: path for path in paths}
    return paths

EXPLANATIONS = __load_explanations(DATA/"explain")

def explain(keyword=None):
    """Print an explanation for an IgSeq concept to stdout"""
    if keyword:
        try:
            path = EXPLANATIONS[keyword]
        except KeyError as err:
            raise IgSeqError(f"key word should be one of {list(EXPLANATIONS.keys())}") from err
        show_text(path)
    else:
        print("Available entries:")
        for key, path in EXPLANATIONS.items():
            with open(path, encoding="UTF8") as f_in:
                preview = next(f_in).rstrip()
            print(f"    {key}  {preview}")
