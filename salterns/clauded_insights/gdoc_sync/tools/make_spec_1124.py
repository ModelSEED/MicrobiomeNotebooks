"""Build gdoc_sync_2026-11-24.json: the 11-20 spec plus the Google Doc changes up to Drive modifiedTime
2026-11-24T16:18:05Z, read through the Docs API (with comments) on 2026-11-24:
  - A. Freiburger's direct edits of 2026-11-23 (22:40Z): Introduction July P24 rewritten again; " in triplicate" added to
    the field-sites paragraph (July P27);
  - A. Freiburger's comment of 2026-11-23 (20:51Z) on the Metagenomics Methods paragraph (July P30);
  - the Docs API suggestion test of 2026-11-24 (suggest.hzsyv3lhfhtb, July P37 "Figure S9" -> "S10"), an open suggestion
    by A. Freiburger; the review's identical P37 edit is marked in_doc_since 2026-11-24 and skipped on this mirror.
Usage: python make_spec_1124.py <read_doc JSON with comments> <out spec>"""
import copy, json, sys
GS = "/home/freiburger/Documents/MicrobiomeNotebooks/salterns/clauded_insights/gdoc_sync/"
old = json.load(open(GS + "gdoc_sync_2026-11-20.json"))
d = json.load(open(sys.argv[1])); api = d["content"] if "content" in d else d
body = api["tabs"][0]["documentTab"]["body"]["content"]
paras = [e["paragraph"] for e in body if "paragraph" in e]
src = open(GS + "tools/make_spec.py").read()
ns = {}
exec(src[src.index("IGNORED ="):src.index("new = copy.deepcopy(old)")], ns)      # segments(), text()
segments, text = ns["segments"], ns["text"]
new = copy.deepcopy(old)
new["source"] = ("Google Doc 'Metabolic Manuscript' (1idf01km9tbuGXXXMU4Hw7BvU0UGpdX2BgShv7BVHcTw), revision "
                 f"{api['revisionId'][:12]}… (Drive modifiedTime 2026-11-24T16:18:05Z), read via the Google Docs API "
                 "(text, suggestions, comments) on 2026-11-24. Changes since the 2026-11-20 spec: A. Freiburger's direct "
                 "edits of 2026-11-23 to July P24 (rewritten) and P27 (' in triplicate'), his comment of 2026-11-23 on "
                 "July P30, and the 2026-11-24 suggestion test (P37).")
# 1) P24 rewritten again (API body paragraph 25); its comment thread AAAB4LFmYcs (Word 12, 13) keeps its anchor
API_P24 = 25
seg24 = segments(paras[API_P24])
t24 = text(seg24)
assert t24.startswith("Which of these hypothesized sources") and "in triplicate" not in t24
assert "such as TMA were demonstrated to fuel methane emissions" in t24
assert "substrate, revealing that restoration efforts should" in t24
th = [c for c in api["comments"] if c["commentId"] == "AAAB4LFmYcs"]
assert th and th[0]["status"] == "OPEN", "the P24 thread changed"
de = new["direct_edits"]
e24 = [e for e in de["edits"] if e["kind"] == "rewrite" and e["p"] == 24]
assert len(e24) == 1
e24[0]["segments"] = seg24
e24[0]["what"] += "; rewritten again 2026-11-23 (its last sentence ends unfinished: '...restoration efforts should ')"
# 2) P27: " in triplicate" (API body paragraph 28); the paragraph also holds Susannah's open "(2021)" -> "(2022)"
p28 = "".join(e["textRun"]["content"] for e in paras[28]["elements"] if "textRun" in e)
assert "122.13340˚W) in triplicate. Water temperature" in p28
de["edits"].append({"p": 27, "kind": "replace", "old": "122.13340˚W). Water", "new": "122.13340˚W) in triplicate. Water",
                    "what": "Field-sites paragraph: ' in triplicate' added after the coordinates, 2026-11-23"})
de["dates"] = de["dates"] + ["2026-11-23"]
de["note"] += (" 2026-11-23: P24's last three sentences are rewritten again (the knowledge gap is now finished: '...such "
               "as TMA were demonstrated to fuel methane emissions in these communities'); P27 gains ' in triplicate'.")
# 3) the suggestion test of 2026-11-24 (open suggestion by A. Freiburger)
sid = "suggest.hzsyv3lhfhtb"
p38 = paras[38]
ins = [e["textRun"]["content"] for e in p38["elements"] if sid in e.get("textRun", {}).get("suggestedInsertionIds", [])]
dele = [e["textRun"]["content"] for e in p38["elements"] if sid in e.get("textRun", {}).get("suggestedDeletionIds", [])]
assert ins == ["10"] and dele == ["9"], (ins, dele)
new["suggestions"]["edits"].append({"p": 37, "kind": "replace", "old": "(see Figure S9)", "new": "(see Figure S10)",
                                    "author": "A. Freiburger", "date": "2026-11-24T16:18:05Z", "suggestion_id": sid,
                                    "note": "Docs API suggestion test (writeMode SUGGEST); the review's own P37 edit"})
new["suggestions"]["note"] += "; plus A. Freiburger's suggestion test of 2026-11-24 (P37)"
# 4) A. Freiburger's comment of 2026-11-23 on the Metagenomics Methods paragraph (July P30). Its Google anchor is the
#    first 603 characters of the paragraph, across citation hyperlinks; the mirror anchors it on the opening words.
c = [x for x in api["comments"] if x["commentId"] == "doco.T54mP-6ngDU4dxnvL2AYgg"]
assert len(c) == 1 and c[0]["status"] == "OPEN" and not c[0].get("replies")
hp = c[0]["headPost"]
assert c[0]["plainTextQuote"].startswith("The 24 metagenomes were quality filtered and trimmed with BBtools")
new["add_comments"].append({"id": "doco.T54mP-6ngDU4dxnvL2AYgg", "anchor": "The 24 metagenomes were quality filtered and "
                            "trimmed with BBtools", "author": "A. Freiburger", "initials": "AF",
                            "date": hp["createTime"][:19] + "Z", "text": hp["content"], "p": 30,
                            "context": "The 24 metagenomes were quality filtered and trimmed with BBtools, assembled",
                            "note": "Google anchor: the first 603 characters of the paragraph (across hyperlinks)"})
new["expected_views"] = {"paragraphs": 492, "note": "unchanged from 2026-11-20"}
json.dump(new, open(sys.argv[2], "w"), indent=1, ensure_ascii=False)
print("P24:", repr(t24[:120]), "...", repr(t24[-80:]))
print("comments:", len(new["add_comments"]), "suggestions:", len(new["suggestions"]["edits"]))
