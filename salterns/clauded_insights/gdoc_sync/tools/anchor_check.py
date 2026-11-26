import ast, json, re, sys
from docx import Document
from docx.oxml.ns import qn
from lxml import etree
CI = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/clauded_insights/"
SAL = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/"
old = Document(SAL + "Metabolic Manuscript - Google Doc mirror (2026-11-03).docx").paragraphs
new = Document(SAL + "Metabolic Manuscript - Google Doc mirror (2026-11-17).docx").paragraphs
def ser(p): return re.sub(r'(<w:(?:ins|del|rPrChange)[^>]*?) w:id="\d+"', r'\1', etree.tostring(p._p).decode())
def direct_text(p):   # what tracked_replace/_locate can see: direct text runs only
    return "".join(r.find(qn("w:t")).text or "" for r in p._p if r.tag == qn("w:r") and r.find(qn("w:t")) is not None)
newser = {}
for j, p in enumerate(new): newser.setdefault(ser(p), []).append(j)
src = open(CI + "make_r232_review_docx.py").read(); tree = ast.parse(src)
main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
idx = set()
for n in ast.walk(main):
    if isinstance(n, ast.Subscript) and isinstance(n.value, ast.Name) and n.value.id == "P" and isinstance(n.slice, ast.Constant):
        idx.add(n.slice.value)
    if isinstance(n, ast.Dict) is False and isinstance(n, ast.DictComp) and ast.unparse(n.key) == "i":
        idx |= set(ast.literal_eval(n.generators[0].iter))
    if isinstance(n, ast.For) and isinstance(n.target, (ast.Name, ast.Tuple)):
        try:
            it = ast.literal_eval(n.iter)
            idx |= {x[0] if isinstance(x, tuple) else x for x in it if isinstance(x if not isinstance(x, tuple) else x[0], int)}
        except Exception: pass
# tracked_math_swap loop items are tuples with str parts; keep ints only
idx = {i for i in idx if isinstance(i, int)}
E = json.load(open(CI + "r232_review_edits_2026-10-26.json"))["edits"]
idx_e = {e["paragraph"] for e in E}
print("indices in main():", sorted(idx)); print("indices in edits JSON:", sorted(idx_e))
print("\n== anchor paragraphs: old mirror P[i] -> position of the identical paragraph in the new mirror")
for i in sorted(idx | idx_e):
    js = newser.get(ser(old[i]), [])
    at_i = "same" if ser(new[i]) == ser(old[i]) else f"now: {new[i].text[:50]!r}"
    print(f"P[{i}]  -> {js if js else 'NO IDENTICAL PARAGRAPH'}   | new P[{i}] {at_i}   | old text {old[i].text[:55]!r}")
print("\n== old strings (tracked_replace / tracked_insert_before in main(); replace/insert_before/comment in edits JSON)")
checks = []
for n in ast.walk(main):
    if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in ("tracked_replace", "tracked_insert_before"):
        k = int(re.search(r"\[(\d+)\]", ast.unparse(n.args[0])).group(1)); s = ast.literal_eval(n.args[1])
        checks.append((f"main L{n.lineno} {n.func.id}", k, s))
for e in E:
    if e["kind"] != "replace_paragraph":
        checks.append((f"edits JSON {e['kind']}", e["paragraph"], e["old"]))
broken = 0
for what, k, s in checks:
    c_old, c_same, c_shift = direct_text(old[k]).count(s), direct_text(new[k]).count(s), direct_text(new[k + 1 if k >= 19 else k]).count(s)
    anywhere = [j for j, p in enumerate(new) if s in direct_text(p)]
    flag = "OK at shifted index" if c_shift == 1 else "BROKEN"
    if c_shift != 1: broken += 1
    if c_shift != 1 or k < 25:
        print(f"{what} P[{k}] {s[:60]!r}: old={c_old} new@same={c_same} new@shifted={c_shift} -> {flag}; found in new at {anywhere}")
print(f"{len(checks)} old strings checked; {broken} no longer match at the shifted index; "
      f"{sum(1 for w,k,s in checks if direct_text(new[k]).count(s)==1)} would still match at the unshifted index")
