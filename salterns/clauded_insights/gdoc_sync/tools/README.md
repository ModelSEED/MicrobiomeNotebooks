# gdoc_sync/tools — scripts used for the 2026-11-17 Google Doc re-sync

Copied from the session scratchpad on 2026-11-19 so the re-sync can be repeated.

| Script | Role |
|---|---|
| `docs_views.py` | Docs API `documents.get` JSON → per-paragraph views, with suggestions rejected (base) and accepted (proposed); produced `../gdoc_views_2026-11-17.json` |
| `diff_views.py` | diff two views files (e.g. 2026-11-03 vs 2026-11-17) → changed / added / removed paragraphs |
| `check_comments.py` | compare the mirror's comment anchors with the Drive export's `<comment_start id=…>` markers |
| `anchor_check.py` | check that the paragraph indices and `tracked_replace` anchors used by `../../make_r232_review_docx.py` still hold on a new mirror |
| `make_spec.py` | build `../gdoc_sync_2026-11-17.json` from the 11-03 spec plus the Docs API read (accepted suggestions, direct edits) |
| `map_comments.py` | map Word comments in a mirror to Google comment threads (Docs API) by author + createTime |
| `mirror_comments.py`, `show.py`, `word_diff.py` | small helpers used during the check (comment listing, paragraph display, word-level diff); no docstrings |

The raw API and Drive dumps were not kept (they were large and session-scoped). Re-read the Doc through the connector tools as `../../../Metabolic Manuscript - R232 review notes.md` ("Google Doc sync 2026-11-17") describes.
