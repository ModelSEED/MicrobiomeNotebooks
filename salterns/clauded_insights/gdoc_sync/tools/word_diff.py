import difflib, re
def wdiff(a, b, ctx=40):
    """compact char-level diff: list of (tag, old_snip, new_snip, left context)"""
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    out = []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal": continue
        out.append(f"   [{tag}] …{a[max(0,i1-ctx):i1]!r}  {a[i1:i2]!r} -> {b[j1:j2]!r}  {a[i2:i2+ctx]!r}…")
    return out
