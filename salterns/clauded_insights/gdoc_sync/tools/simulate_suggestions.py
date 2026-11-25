"""Dry run of review_to_suggestions: apply a plan's operations (SUGGEST semantics: deletions mark, insertions add) to
the Docs API read and compare each touched paragraph's proposed view (all suggestions accepted) with the review's
accepted view. Usage: python simulate_suggestions.py <review.docx> <read_doc JSON> <plan.json>"""
import json, sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import review_to_suggestions as R
from docx import Document

review, apip, planp = sys.argv[1:4]
api = json.load(open(apip)); api = api.get("content", api)
P = json.load(open(planp))
# document model: list of [index, char, state] over the whole body text runs (state: base / ins / del)
chars = []
AP = R.api_paragraphs(api)
for a in AP:
    for e in a["paragraph"]["elements"]:
        if "textRun" not in e:
            chars.append([e["startIndex"], "\x00", "obj"]); continue
        tr = e["textRun"]; i = e["startIndex"]
        st = "ins" if tr.get("suggestedInsertionIds") else ("del" if tr.get("suggestedDeletionIds") else "base")
        for ch in tr["content"]:
            chars.append([i, ch, st]); i += R.u16(ch)
_, reqs = R.chunk_requests(P, 0, 10 ** 9)
for rq in reqs:
    if "deleteContentRange" in rq:
        r = rq["deleteContentRange"]["range"]
        for c in chars:
            if r["startIndex"] <= c[0] < r["endIndex"] and c[2] != "ins":
                c[2] = "del"
    elif "insertText" in rq:
        at, t = rq["insertText"]["location"]["index"], rq["insertText"]["text"]
        pos = next(n for n, c in enumerate(chars) if c[0] >= at)
        for c in chars[pos:]:
            c[0] += len(t)
        chars[pos:pos] = [[at + n, ch, "ins"] for n, ch in enumerate(t)]
prop = "".join(c[1] for c in chars if c[2] in ("base", "ins") and c[1] != "\x00")
doc = Document(review)
acc = []
for p in doc.element.body.iter(f"{R.W}p"):
    tk, _ = R.tokens(p)
    acc.append("".join(ch for ch, st, f in tk if st in ("base", "cins", "oins")))
# compare as paragraph lists split on newlines (review paragraphs that are pictures/tables are not in the Doc)
docp = [x for x in prop.split("\n")]
import difflib
sm = difflib.SequenceMatcher(None, [x.strip() for x in docp if x.strip()], [x.strip() for x in acc if x.strip()], autojunk=False)
bad = [(t, i2 - i1, j2 - j1) for t, i1, i2, j1, j2 in sm.get_opcodes() if t != "equal"]
print("paragraph blocks differing:", len(bad))
for t, i1, i2, j1, j2 in sm.get_opcodes():
    if t == "equal":
        continue
    a = [x.strip() for x in docp if x.strip()][i1:i2]; b = [x.strip() for x in acc if x.strip()][j1:j2]
    if any(x[:12] in ("GTDB-Tk Fami", "0.205579") for x in a) or (not a and len(b) > 2 and len(b[0]) < 40):
        print("== table block", t, len(a), len(b)); continue
    for x, y in zip(a + [""] * (len(b) - len(a)), b + [""] * (len(a) - len(b))):
        d = [(tt, repr(x[i1:i2]), repr(y[j1:j2])) for tt, i1, i2, j1, j2 in difflib.SequenceMatcher(None, x, y, autojunk=False).get_opcodes() if tt != "equal"]
        print("==", t, "|", y[:70] or x[:70], "|", d[:4])
