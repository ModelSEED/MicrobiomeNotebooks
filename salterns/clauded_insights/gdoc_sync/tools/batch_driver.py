"""Drive the suggestion batch in chunks, from the end of the Doc backward (review_to_suggestions.py makes the plan).

  chunks  <plan.json> <out_dir> [--first N] [--max N]   split the operations into chunk_NN.json (+ chunk_NN.req.json,
                                                        the compact request list to send)
  precheck <initial read> <fresh read> <chunk.json>     the Doc below the chunk's highest index is unchanged since
                                                        the plan's read; prints the fresh revisionId
  verify  <initial read> <fresh read> <out_dir> <n>     the Doc equals the initial read plus chunks 1..n, in both the
                                                        base view and the proposed view (text), and the inserted
                                                        ranges carry the planned italics
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import review_to_suggestions as R


def load(p):
    d = json.load(open(p))
    return d.get("content", d)


def chars_of(api):
    out = []
    for a in R.api_paragraphs(api):
        for e in a["paragraph"]["elements"]:
            if "textRun" not in e:
                out.append([e["startIndex"], "\x00", "obj", {}]); continue
            tr = e["textRun"]; i = e["startIndex"]
            st = "ins" if tr.get("suggestedInsertionIds") else ("del" if tr.get("suggestedDeletionIds") else "base")
            ts = dict(tr.get("textStyle", {}))
            for v in tr.get("suggestedTextStyleChanges", {}).values():     # the style after the suggestions
                if v.get("textStyleSuggestionState", {}).get("italicSuggested"):
                    ts["italic"] = v["textStyle"].get("italic", False)
            for ch in tr["content"]:
                out.append([i, ch, st, ts]); i += R.u16(ch)
    return out


def chunks(plan_p, out, first=4000, mx=9000, start=1):
    P = json.load(open(plan_p))
    sel, _ = R.chunk_requests(P, 0, 10 ** 9)
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    groups, cur, size = [], [], 0
    limit = first
    for o in sel:
        n = len(json.dumps(R.op_requests(o), ensure_ascii=False, separators=(",", ":")))
        if cur and size + n > limit:
            groups.append(cur); cur, size = [], 0; limit = mx
        cur.append(o); size += n
    if cur:
        groups.append(cur)
    for k, g in enumerate(groups, start):
        reqs = [r for o in g for r in R.op_requests(o)]
        meta = {"chunk": k, "max_index": max(o["index"] for o in g), "min_index": min(o["index"] for o in g),
                "sections": sorted({o["section"] for o in g}), "n_ops": len(g), "n_requests": len(reqs), "ops": g}
        json.dump(meta, open(out / f"chunk_{k:02d}.json", "w"), indent=1, ensure_ascii=False)
        json.dump(reqs, open(out / f"chunk_{k:02d}.req.json", "w"), ensure_ascii=False, separators=(",", ":"))
        print(k, meta["max_index"], meta["min_index"], len(g), len(reqs),
              len(json.dumps(reqs, ensure_ascii=False, separators=(",", ":"))), meta["sections"][:3])


def precheck(init_p, fresh_p, chunk_p):
    a, b = load(init_p), load(fresh_p)
    X = json.load(open(chunk_p))["max_index"]
    ca = [c[:3] for c in chars_of(a) if c[0] < X]
    cb = [c[:3] for c in chars_of(b) if c[0] < X]
    assert ca == cb, "the Doc changed below the chunk; re-plan"
    print("unchanged up to index", X, "| revisionId", b["revisionId"])


def apply(chars, reqs):
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
            chars[pos:pos] = [[at + n, ch, "ins", {}] for n, ch in enumerate(t)]
        elif "updateTextStyle" in rq:
            r = rq["updateTextStyle"]["range"]; ts = rq["updateTextStyle"]["textStyle"]
            for c in chars:
                if r["startIndex"] <= c[0] < r["endIndex"]:
                    c[3] = dict(c[3], **{"italic": ts["italic"]} if "italic" in ts else {})
    return chars


def verify(init_p, fresh_p, out, n):
    """out: a directory of chunk_NN.req.json (chunks 1..n), or a text file listing the sent request files in order."""
    a, b = load(init_p), load(fresh_p)
    exp = chars_of(a)
    files = ([Path(out) / f"chunk_{k:02d}.req.json" for k in range(1, n + 1)] if Path(out).is_dir()
             else [Path(x) for x in open(out).read().split()][:n])
    for fp in files:
        exp = apply(exp, json.load(open(fp)))
    got = chars_of(b)
    views = lambda cs, keep: "".join(c[1] for c in cs if c[2] in keep and c[1] != "\x00")
    ok = True
    for name, keep in (("base", ("base", "del")), ("proposed", ("base", "ins"))):
        e, g = views(exp, keep), views(got, keep)
        if e != g:
            ok = False
            i = next(j for j in range(min(len(e), len(g))) if e[j] != g[j]) if any(
                e[j] != g[j] for j in range(min(len(e), len(g)))) else min(len(e), len(g))
            print(f"{name} view differs at char {i}: expected {e[i - 60:i + 60]!r}\n  got {g[i - 60:i + 60]!r}")
    # italics of the inserted characters (by order of the inserted text)
    ei = [bool(c[3].get("italic")) for c in exp if c[2] == "ins"]
    gi = [bool(c[3].get("italic")) for c in got if c[2] == "ins"]
    et = "".join(c[1] for c in exp if c[2] == "ins"); gt = "".join(c[1] for c in got if c[2] == "ins")
    if et == gt and ei != gi:
        bad = sum(x != y for x, y in zip(ei, gi))
        print(f"italics differ on {bad} of {len(ei)} inserted characters"); ok = False
    print("OK" if ok else "MISMATCH", "| suggestions in the Doc now:", sum(1 for c in got if c[2] == "ins"), "inserted chars")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "chunks":
        args = sys.argv[2:]
        kw = {}
        if "--first" in args:
            kw["first"] = int(args[args.index("--first") + 1])
        if "--max" in args:
            kw["mx"] = int(args[args.index("--max") + 1])
        if "--start" in args:
            kw["start"] = int(args[args.index("--start") + 1])
        chunks(args[0], args[1], **kw)
    elif cmd == "precheck":
        precheck(*sys.argv[2:5])
    elif cmd == "verify":
        verify(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
