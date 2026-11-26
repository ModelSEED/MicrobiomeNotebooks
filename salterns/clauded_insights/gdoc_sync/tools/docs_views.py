"""Docs API documents.get JSON -> per-paragraph views (suggestions rejected = base, accepted = proposed)."""
import json, sys, collections

def load(path):
    d = json.load(open(path))
    return d["content"] if "content" in d else d

def elem_text(e):
    if "textRun" in e: return e["textRun"]["content"], e["textRun"]
    if "inlineObjectElement" in e: return "[IMG]", e["inlineObjectElement"]
    if "equation" in e: return "[EQ]", e["equation"]
    k = [x for x in e if x not in ("startIndex", "endIndex")][0]
    return None, (k, e[k])

def para_views(p, unknown):
    base, prop = [], []
    for e in p["elements"]:
        t, obj = elem_text(e)
        if t is None:
            unknown[obj[0]] += 1
            continue
        ins = bool(obj.get("suggestedInsertionIds")); dele = bool(obj.get("suggestedDeletionIds"))
        if not ins: base.append(t)
        if not dele: prop.append(t)
    return "".join(base), "".join(prop)

def views(c):
    body = c["tabs"][0]["documentTab"]["body"]["content"]
    out, unknown = [], collections.Counter()
    for el in body:
        if "paragraph" in el:
            b, p = para_views(el["paragraph"], unknown)
            out.append({"i": len(out), "where": "body", "base": b, "proposed": p, "_start": el.get("startIndex")})
        elif "table" in el:
            for ri, row in enumerate(el["table"]["tableRows"]):
                for ci, cell in enumerate(row["tableCells"]):
                    for ce in cell["content"]:
                        if "paragraph" in ce:
                            b, p = para_views(ce["paragraph"], unknown)
                            out.append({"i": len(out), "where": f"table r{ri}c{ci}", "base": b, "proposed": p, "_start": ce.get("startIndex")})
                        else:
                            unknown["cell:" + [k for k in ce if k not in ('startIndex','endIndex')][0]] += 1
    return out, unknown

if __name__ == "__main__":
    c = load(sys.argv[1])
    v, unk = views(c)
    print(len(v), dict(unk), file=sys.stderr)
    json.dump(v, open(sys.argv[2], "w"), indent=1, ensure_ascii=False)
