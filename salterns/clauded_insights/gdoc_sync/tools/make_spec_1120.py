"""Build gdoc_sync_2026-11-20.json: the 11-17 spec plus A. Freiburger's direct edits of 2026-11-20 (Drive modifiedTime
2026-11-20T13:56:07Z) to three Introduction paragraphs (July P21, P23, P24), read through the Docs API on 2026-11-20.
Usage: python make_spec_1120.py <read_doc JSON> <out spec>"""
import copy, json, sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
GS = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/clauded_insights/gdoc_sync/"
old = json.load(open(GS + "gdoc_sync_2026-11-17.json"))
d = json.load(open(sys.argv[1])); api = d["content"] if "content" in d else d
body = api["tabs"][0]["documentTab"]["body"]["content"]
paras = [e["paragraph"] for e in body if "paragraph" in e]
src = open(GS + "tools/make_spec.py").read()
ns = {}
exec(src[src.index("IGNORED ="):src.index("new = copy.deepcopy(old)")], ns)      # segments(), text()
segments, text = ns["segments"], ns["text"]
new = copy.deepcopy(old)
new["source"] = ("Google Doc 'Metabolic Manuscript' (1idf01km9tbuGXXXMU4Hw7BvU0UGpdX2BgShv7BVHcTw), revision "
                 f"{api['revisionId'][:12]}… (Drive modifiedTime 2026-11-20T13:56:07Z), read via the Google Docs API on "
                 "2026-11-20. Changes since the 2026-11-17 spec: A. Freiburger's direct edits of 2026-11-20 to July P21, P23 "
                 "and P24 (Introduction); no suggestion or comment changes in the API read.")
API_OF = {21: 23, 23: 24, 24: 25}                     # July index -> API body-paragraph index (after the 11-17 move)
seg = {p: segments(paras[a]) for p, a in API_OF.items()}
assert text(seg[21]).startswith("\tThe biogeochemistry metabolism") and "(38)" in text(seg[21])
assert text(seg[23]).startswith("Previous studies") and "These salt-adapted methanogens are" in text(seg[23])
assert text(seg[24]).startswith("Which of these hypothesized sources")
de = new["direct_edits"]
for e in de["edits"]:
    if e["kind"] == "rewrite" and e["p"] in (21, 23):
        e["segments"] = seg[e["p"]]
        e["what"] += "; rewritten again 2026-11-20"
de["edits"].insert(next(i for i, e in enumerate(de["edits"]) if e["kind"] != "rewrite"),
                   {"p": 24, "kind": "rewrite", "what": "Knowledge-gap paragraph ('Which of these hypothesized sources...'), "
                    "rewritten 2026-11-20 (its second sentence is unfinished in the Google Doc: '...such as TMA were .')",
                    "segments": seg[24]})
de["dates"] = de["dates"] + ["2026-11-20"]
de["note"] += (" 2026-11-20: P21 gains the citation (38) after 'hydrogenotrophic methanogenesis'; P23's last three "
               "sentences are rewritten (the '(38)' citation moves to P21); P24's first two sentences are rewritten.")
for e in de["edits"]:
    if e["kind"] == "rewrite" and e["p"] == 23:
        e.pop("reanchor_comments", None)                 # Word id 7 ("methanogen" thread) resolved 2026-11-20
    if e["kind"] == "rewrite" and e["p"] == 24:
        e["reanchor_comments"] = [{"word_ids": ["12", "13"], "thread": "AAAB4LFmYcs",
                                   "anchor": "revealing that restoration efforts should",
                                   "context": "substrate, revealing that restoration efforts should"}]
rc = new["remove_comments"]
rc["word_ids"] = rc["word_ids"] + ["7", "9", "10", "11"]
new["add_replies"] = [r for r in new["add_replies"] if r["parent_word_id"] != "10"]   # Susannah's 11-03 reply (Word id 84)
rc["note"] += (" 2026-11-20: resolved by A. Freiburger (Docs API commentAction RESOLVE): Word id 7, 'Although \"methanogen\" "
               "technically just means...' (AAAB50XBNtQ, 03:55Z); Word id 9, 'Most abundant genes in this environment "
               "overall?...' (AAAB50XBNug, 03:55Z); Word ids 10, 11, 84, 'what was the evidence for betaine being the "
               "methanogenic substrate?' with Clifton's and Susannah's replies (AAAB4LFmYcE, 13:54Z; its anchor sentence was "
               "deleted in the P24 rewrite).")
new["expected_views"] = {"paragraphs": 492, "note": "unchanged from 2026-11-17"}
json.dump(new, open(sys.argv[2], "w"), indent=1, ensure_ascii=False)
for p, v in seg.items(): print(p, len(v), repr(text(v)[:90]))
