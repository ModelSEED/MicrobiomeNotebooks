#!/bin/bash
# Rebuild the salterns 2026 combined ModelSEED template set.
#
# This is the exact sequence that produced the four JSON files in modelling_2026/templates/
# on 2026-11-20 (branch salterns-2026-combined, head 9696e5f in ~/Documents/ModelSEEDTemplates).
# A rerun gives byte-identical JSON (same md5s as PROVENANCE.md) but new commit hashes.
#
# Inputs, all read-only:
#   ~/Documents/ModelSEEDTemplates               codif-curated (b6fe86b, PR #14), 7fd1ab1 (PR #15),
#                                                0f6a5f6 (PR #16)
#   ~/Documents/ModelSEEDTemplates-hppase        e57ea75 (PR #17; this commit is not in the main clone)
# The three apply scripts are copied verbatim from those commits into build/source_scripts/.
#
# Usage: build_combined.sh [WORKTREE_DIR] [BRANCH] [OUT_DIR]
set -euo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
REPO=/home/freiburger/Documents/ModelSEEDTemplates
HPPASE=/home/freiburger/Documents/ModelSEEDTemplates-hppase
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${3:-$(dirname "$HERE")}"                 # default: modelling_2026/templates
WT="${1:-/home/freiburger/Documents/ModelSEEDTemplates-salterns2026}"
BR="${2:-salterns-2026-combined}"
S="$HERE/source_scripts"
L="${LOG_DIR:-$HERE/logs}"
mkdir -p "$S" "$L"

# 0. the apply scripts, verbatim from the PR commits
git -C "$REPO"   show 7fd1ab1:validation/apply_draft_loop_fix.py        > "$S/apply_draft_loop_fix.py"
git -C "$REPO"   show 0f6a5f6:validation/apply_thiolase_direction.py    > "$S/apply_thiolase_direction.py"
if [ -d "$HPPASE" ]; then   # e57ea75 exists locally only in this clone; fall back to the kept copy
  git -C "$HPPASE" show e57ea75:validation/apply_proton_pyrophosphatase.py > "$S/apply_proton_pyrophosphatase.py"
fi
# the copies must be the ones recorded in PROVENANCE.md
(cd "$S" && md5sum -c - <<'MD5'
6d59d7dcf87b918e4e8aabfb9fb85842  apply_draft_loop_fix.py
72d204377938dcb5822bed40463e5dbe  apply_thiolase_direction.py
981fde3994d6aad357d3bbea7e56f9e8  apply_proton_pyrophosphatase.py
MD5
)

# 1. new branch + worktree on codif-curated (b6fe86b = PR #14 head; carries #9, #11, #12, #13)
test "$(git -C "$REPO" rev-parse codif-curated)" = "b6fe86bd880ef44ccb2f1104a2b6f4c163dd74df"
git -C "$REPO" worktree add -b "$BR" "$WT" b6fe86bd880ef44ccb2f1104a2b6f4c163dd74df
cd "$WT"

# 2. PR #15 (fix-draft-atp-loops, 7fd1ab1): seven draft-loop reactions forward
$PY "$S/apply_draft_loop_fix.py" --templates "$WT/templates/v7.0" > "$L/01_apply_draft_loop_fix.json"
$PY "$S/apply_draft_loop_fix.py" --templates "$WT/templates/v7.0" --check > "$L/01_check.json"
git add templates/v7.0
git commit -q -m "Replay #15 (fix-draft-atp-loops, 7fd1ab1) on codif-curated: seven draft-loop reactions forward"

# 3. PR #16 (atp-cycle-directions, 0f6a5f6): rxn00598 thiolysis '<'
$PY "$S/apply_thiolase_direction.py" --templates "$WT/templates" > "$L/02_apply_thiolase_direction.json"
$PY "$S/apply_thiolase_direction.py" --templates "$WT/templates" --check > "$L/02_check.json"
git add templates
git commit -q -m "Replay #16 (atp-cycle-directions, 0f6a5f6): rxn00598 3-oxoadipyl-CoA thiolase to thiolysis '<'"

# 4. PR #17 (hppase-pyrophosphate-proton-pump, e57ea75): rxn42569 + tcpx01431
$PY "$S/apply_proton_pyrophosphatase.py" --templates "$WT/templates" > "$L/03_apply_proton_pyrophosphatase.json"
$PY "$S/apply_proton_pyrophosphatase.py" --templates "$WT/templates" --check > "$L/03_check.json"
git add templates
git commit -q -m "Replay #17 (hppase-pyrophosphate-proton-pump, e57ea75): rxn42569 H+-PPase for orphan role ftr05360"

# 5. export the four v7.0 files from the committed tree
for f in GramNegModelTemplateV7.json GramPosModelTemplateV7.json ArchaeaTemplateV6.json Core-V6.json; do
  git show "$BR:templates/v7.0/$f" > "$OUT/$f"
done
(cd "$OUT" && md5sum GramNegModelTemplateV7.json GramPosModelTemplateV7.json ArchaeaTemplateV6.json Core-V6.json)
