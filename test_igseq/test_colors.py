from igseq import colors
from .util import TestBase


class TestColors(TestBase):
    """Tests for color-related helper functions."""

    TRIOS = [
        ([255,   0,   0], "#ff0000"),
        ([0,     0,   0], "#000000"),
        ([0,   136,   0], "#008800"),
        ([255, 255, 255], "#ffffff"),
        ]

    TEXTS = [
        ("#ff0000", [255, 0, 0]),
        ("#FF0000", [255, 0, 0]),
        ("FF0000",  [255, 0, 0]),
        ("#008800", [0, 136, 0]),
        ("#f00",    [255, 0, 0]),
        ("#F00",    [255, 0, 0]),
        ("#080",    [0, 136, 0]),
        ("f00",     [255, 0, 0]),
        ]

    SCALES = [
        # Two colors averaged, no scaling
        (([[255, 0, 0], [0, 0, 255]], 0), [127, 0, 127]),
        # Two colors averaged, of 2 total, scales to black
        (([[255, 0, 0], [0, 0, 255]], 2), [0, 0, 0]),
        # one color of two, stays the same
        (([[255, 0, 0]], 2), [255, 0, 0]),
        # no colors = black by definition
        (([], 0), [0, 0, 0]),
        (([], 2), [0, 0, 0]),
        # two colors of three, averaged + scaled
        (([[255, 0, 0], [0, 0, 255]], 3), [91, 0, 91]),
        ]

    def test_merge_colors(self):
        """Test blending colors together."""
        for case in self.__class__.SCALES:
            with self.subTest(case=case):
                self.assertEqual(colors.merge_colors(case[0][0], case[0][1]), case[1])

    def test_color_str_to_trio(self):
        """Test converting color text codes to integer trios."""
        for case in self.__class__.TEXTS:
            with self.subTest(case=case):
                self.assertEqual(colors.color_str_to_trio(case[0]), case[1])

    def test_color_trio_to_str(self):
        """Test converting integer trios to color text codes."""
        for case in self.__class__.TRIOS:
            with self.subTest(case=case):
                self.assertEqual(colors.color_trio_to_str(case[0]), case[1])
