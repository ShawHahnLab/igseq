
# from SONAR's display_tree, which in turn got these from the iwanthue tool
COLORS = ["#BE4229", "#E74721", "#8C431C", "#CB6A27", "#E98C25", "#946E13", "#D0A620", "#8B8A22", "#A9B71D", "#555C10", "#839A21", "#82B531", "#4D8426", "#55C42A", "#34992A", "#2B6B2E", "#4FC456", "#33B66D", "#4296CB", "#5A8CE2", "#3E5988", "#656CE2", "#524EA0", "#8F83CC", "#A57CE4", "#8E46AD", "#C056EB", "#CA6BE4", "#7B4D87", "#D186D7"]

def merge_colors(colors, scale=0):
    """Take an average of a list of colors and shift toward black.

    More colors results in a darker result, up to the integer value given for
    scale.  If scale is less than the number of colors this scaling is skipped.
    """
    result = [0, 0, 0]
    if not colors:
        return result
    if len(colors) == 1:
        return colors[0]
    for color in colors:
        for idx in range(3):
            result[idx] += color[idx]
    # not quite right, should rotate, really, not move directly toward the
    # middle... but it'll do for now
    if scale < len(colors):
        scaling = 1
    else:
        scaling = ((scale - len(colors))/scale)**0.3
    for idx in range(3):
        result[idx] = result[idx] / len(colors)
        result[idx] = int(result[idx] * scaling)
    return result

def color_str_to_trio(color_txt):
    """Convert hex color string to trio of 0:255 ints."""
    color_txt = color_txt.removeprefix("#")
    # e.g. "ff0000"
    if len(color_txt) == 6:
        color = [int(color_txt[idx:(idx+2)], 16) for idx in range(0, 6, 2)]
    # e.g. "f00" = "ff0000"
    elif len(color_txt) == 3:
        color = [int(color_txt[idx:(idx+1)]*2, 16) for idx in range(0, 3)]
    else:
        raise ValueError
    return color

def color_trio_to_str(color):
    """Convert trio of 0:255 ints to hex color string."""
    return "#" + "".join([f"{c:02x}" for c in color])
