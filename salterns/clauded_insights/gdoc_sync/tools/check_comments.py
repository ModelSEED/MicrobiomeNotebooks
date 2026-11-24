"""Every Word comment in the mirror must be an OPEN Google thread with export markers, and its anchored text must equal
the export's marker-enclosed text (markdown stripped, whitespace ignored), in the base or the accepted view."""
import sys, json, re
sys.path.insert(0, sys.argv[4])
from map_comments import api_posts, export_anchors, mdnorm, mirror_comments
mirror, api, export = sys.argv[1:4]
comments, posts = api_posts(api)
status = {t["commentId"]: t["status"] for t in comments}
anchors = export_anchors(export)
doc, meta, res = mirror_comments(mirror)
ok = bad = 0; threads = set()
def wn(s): return re.sub(r"\s+", "", s.replace(" ", " "))
for wid, (author, date, text) in meta.items():
    hits = posts.get((author, date), [])
    assert len(hits) == 1, (wid, author, date, hits)
    tid, kind, st, quote = hits[0]
    threads.add(tid)
    segs = anchors.get(tid)
    r = res.get(wid)
    if st != "OPEN" or segs is None or not r:
        print("PROBLEM", wid, tid, st, segs is not None, bool(r)); bad += 1; continue
    exp = mdnorm("".join(segs))
    mb = wn("".join(x["base"] for x in r)); ma = wn("".join(x["acc"] for x in r))
    if exp in (mb, ma):
        ok += 1
    else:
        print("MISMATCH", wid, tid, "export:", exp[:80], "| mirror base:", mb[:80], "| acc:", ma[:80]); bad += 1
open_anch = {t for t, s in status.items() if s == "OPEN" and t in anchors}
print(f"{ok} comments match the export anchors, {bad} problems; threads in mirror {len(threads)}, open anchored threads {len(open_anch)}, missing {open_anch - threads}")
