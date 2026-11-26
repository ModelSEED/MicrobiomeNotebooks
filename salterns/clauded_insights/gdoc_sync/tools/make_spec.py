"""Build gdoc_sync_2026-11-17.json from the 11-03 spec + the Docs API read of 2026-11-19 (revision of 2026-11-17)."""
import copy, json, sys
GS = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/clauded_insights/gdoc_sync/"
old = json.load(open(GS + "gdoc_sync_2026-11-03.json"))
api = json.load(open("raw/read_doc_comments_2026-11-19.json"))["content"]
body = api["tabs"][0]["documentTab"]["body"]["content"]
paras = [e["paragraph"] for e in body if "paragraph" in e]          # body paragraphs (tables excluded), API order

IGNORED = {"weightedFontFamily", "fontSize", "foregroundColor", "backgroundColor"}
def segments(p):
    """API paragraph -> list of segments: plain string, or {"t", "i"/"b"/"sup"/"sub"/"link"}."""
    out = []
    for e in p["elements"]:
        tr = e["textRun"]; ts = tr.get("textStyle", {})
        assert not tr.get("suggestedInsertionIds") and not tr.get("suggestedDeletionIds"), tr
        extra = set(ts) - IGNORED - {"italic", "bold", "baselineOffset", "link", "underline"}
        assert not extra, extra
        st = {}
        if ts.get("italic"): st["i"] = True
        if ts.get("bold"): st["b"] = True
        if ts.get("baselineOffset") == "SUPERSCRIPT": st["sup"] = True
        if ts.get("baselineOffset") == "SUBSCRIPT": st["sub"] = True
        if "link" in ts: st["link"] = ts["link"]["url"]
        t = tr["content"].rstrip("\n") if e is p["elements"][-1] else tr["content"]
        if not t: continue
        if out and (out[-1] if isinstance(out[-1], dict) else {"t": out[-1]}).keys() - {"t"} == st.keys() and \
           all((out[-1] if isinstance(out[-1], dict) else {}).get(k) == v for k, v in st.items()) and "link" not in st:
            if isinstance(out[-1], dict): out[-1]["t"] += t
            else: out[-1] += t
            continue
        out.append(dict(t=t, **st) if st else t)
    return out

def text(segs): return "".join(s if isinstance(s, str) else s["t"] for s in segs)

new = copy.deepcopy(old)
new["source"] = ("Google Doc 'Metabolic Manuscript' (1idf01km9tbuGXXXMU4Hw7BvU0UGpdX2BgShv7BVHcTw), revision "
                 f"{api['revisionId'][:12]}… (Drive modifiedTime 2026-11-17T23:20:24Z), read via the Google Docs API (text, "
                 "suggestions, comments) and the Drive API (export with comment markers) on 2026-11-19. Changes since the "
                 "2026-11-03 spec were made by A. Freiburger on 2026-11-16 and 2026-11-17; the 2026-11-05 comment threads are "
                 "unchanged.")
# 1) suggestions accepted since 11-03: move out of the open suggestions
acc_ids = {8: ("suggest.xu2p6datvpgo", "2026-11-16T19:03:54Z"), 23: ("suggest.1yr1msb8n3kp", "2026-11-16T19:04:02Z"),
           24: ("suggest.f3bf6e9xl761", "2026-11-17T23:04:00Z")}
keep, accepted = [], []
for e in old["suggestions"]["edits"]:
    if e["p"] in acc_ids:
        a = dict(e); a["suggestion_id"], a["accepted"] = acc_ids[e["p"]]; accepted.append(a)
    else:
        keep.append(e)
assert len(accepted) == 3 and len(keep) == 19
new["suggestions"]["edits"] = keep
new["suggestions"]["note"] = ("Susannah Tringe's suggestions of 2026-11-03 that are still open on 2026-11-17 (19 of the 22 "
                              "edits; the other 3 were accepted, see accepted_suggestions)")
new["accept_all_revisions"] = {
    "note": ("Susannah Tringe's suggestions of 2026-05-07 in these paragraphs (tracked changes in the 2026-07-02 export) were "
             "all accepted by A. Freiburger on 2026-11-17 (24 Google suggestion ids, below). The only comment inside them "
             "(Word id 8, on the deleted 'minerals') loses its anchor; see remove_comments."),
    "paragraphs": {"23": {"author": "Susannah Tringe", "elements": 23}, "24": {"author": "Susannah Tringe", "elements": 12}},
    "suggestion_ids": ["suggest.h61t15g217q5", "suggest.9bfb3ahikf2m", "suggest.h1e0om5l46m5", "suggest.fjnmu89awsau",
                       "suggest.50412fdzcn4n", "suggest.i0ugkacc87ya", "suggest.1j3nqri4m76k", "suggest.q70zgnlnqvbb",
                       "suggest.pv7qmpyqe1wz", "suggest.1nefknp2byag", "suggest.dz362mq6mbfm", "suggest.e6gmaq8fuhms",
                       "suggest.nkramn86chki", "suggest.rrqb27u47br5", "suggest.o1lkvctbw9jf", "suggest.lo4gg2882sqb",
                       "suggest.3tw7nlldn3pn", "suggest.1sgmzl9o8cc4", "suggest.5l20hh3bvjua", "suggest.n09vpoejguhm",
                       "suggest.yzc3hb6t8w1d", "suggest.i6rs7n7dyfq", "suggest.vw4cjtrrmkgn", "suggest.71kp8wx47m6a"]}
new["accepted_suggestions"] = {
    "note": "3 of Susannah Tringe's 2026-11-03 suggestions, accepted by A. Freiburger; applied as plain (untracked) text",
    "edits": accepted}
# 2) comments
new["remove_comments"] = {
    "note": ("Word ids 4, 5: thread resolved by A. Freiburger on 2026-08-27 ('does this mean a line on log-scale axes?' + "
             "Clifton's reply). Word ids 2, 3: trehalose thread (AAAB2s-yGXQ) resolved by A. Freiburger on 2026-11-17. Word id "
             "6: 'Need to confirm that this reference supports my statement about sulfate' (AAAB50XBNtM, on '15') resolved "
             "by A. Freiburger on 2026-11-17. Word id 8: 'Methylphosphonates are also organic molecules, not minerals' "
             "(AAAB50XBNuI) is still open but its anchor 'minerals' was deleted when the suggestion was accepted on "
             "2026-11-17, so Google's exports drop it (12 open threads are now unanchored)."),
    "word_ids": ["4", "5", "2", "3", "6", "8"]}
# 3) direct edits (no suggestion record): abstract split + rewrite, Introduction restructure, P23 head rewrite
P = {k: segments(paras[k]) for k in (18, 19, 21, 22, 23, 24)}
new["direct_edits"] = {
    "note": ("Direct edits in the Google Doc on 2026-11-16/17 (A. Freiburger's session; the Docs API keeps no suggestion record "
             "for them, except his own abstract suggestions suggest.f6bovdvlwphd and suggest.99o1bf1fjf34, accepted "
             "2026-11-16, and suggest.gkv4ukhhru9, rejected). Paragraph indices are those of the July export. Each rewrite "
             "replaces the paragraph's runs with the Docs API runs (italic, sub/superscript, Zotero citation links); "
             "'insert_after' and 'move_before' run last, after all index-based steps."),
    "author": "A. Freiburger", "dates": ["2026-11-16", "2026-11-17"],
    "edits": [
        {"p": 18, "kind": "rewrite", "what": "Abstract, first part (ends at '...underlying biogeochemistry.')", "segments": P[18]},
        {"p": 18, "kind": "insert_after", "what": "Abstract, second part (new paragraph, same paragraph format)", "segments": P[19]},
        {"p": 20, "kind": "rewrite", "what": "Introduction paragraph 1 (salt habitats, salt-in/salt-out; absorbs '> 7% salt (15)' and the compatible-solute list)", "segments": P[21]},
        {"p": 22, "kind": "rewrite", "what": "Compatible solutes (GB, trehalose); now the 2nd Introduction paragraph", "segments": P[22]},
        {"p": 21, "kind": "rewrite", "what": "'The biogeochemistry metabolism of organisms...' (shortened); now the 3rd Introduction paragraph", "segments": P[23]},
        {"p": 22, "kind": "move_before", "before": 21, "what": "the compatible-solute paragraph moves ahead of the biogeochemistry paragraph"},
        {"p": 23, "kind": "rewrite", "what": "Methane-salinity paragraph: Susannah's suggestions accepted, then its first two sentences rewritten",
         "segments": P[24],
         "reanchor_comments": [{"word_ids": ["7"], "thread": "AAAB50XBNtQ", "anchor": "methanogens",
                                "context": "These salt-adapted methanogens are"}]},
    ]}
new["expected_views"] = {"paragraphs": 492, "note": "the abstract split adds one paragraph (491 in the 11-03 views)"}
json.dump(new, open(sys.argv[1], "w"), indent=1, ensure_ascii=False)
for k, v in P.items(): print(k, len(v), repr(text(v)[:70]))
