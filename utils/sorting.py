"""
Natural sorting of the lists of strings containing floats (e.g. 'XX_1.0_YY', 'XX_1.05_YY')
Implementation taken from https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
"""


def atof(text):
    try:
        retval = float(text)
    except ValueError:
        retval = text
    return retval


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    float regex comes from https://stackoverflow.com/a/12643073/190597
    """
    return [atof(c) for c in re.split(r"[+-]?([0-9]+(?:[.][0-9]*)?|[.][0-9]+)", text)]
