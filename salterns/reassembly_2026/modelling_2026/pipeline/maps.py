"""Annotation-term -> ModelSEED reaction maps for the pre-registered annotation comparison (arms B-D).

Re-runnable. Reads only:
  - kegg_cache/link_reaction_ko_<date>.tsv  (KEGG REST https://rest.kegg.jp/link/reaction/ko, cached with its date)
  - ~/Documents/ModelSEEDDatabase Biochemistry/Aliases (KEGG and BiGG reaction aliases, Enzyme Class table),
    read at the checkout's current commit (recorded; never switched)
  - the four combined templates (template reaction ids and role names), to report how much of each map lands in them

Writes maps/term_maps.json:
  ko2rxn   KO -> [ModelSEED rxn ids]   (KO -> KEGG R (KEGG link) -> rxn (ModelSEED 'KEGG' aliases))
  ec2rxn   EC -> [rxn ids]             (ModelSEED 'Enzyme Class' table; complete 4-level ECs only)
  bigg2rxn BiGG reaction id -> [rxn]   (ModelSEED 'BiGG' aliases)
plus provenance (file md5s, DB commit, dates) and coverage against each template.
usage: python maps.py
"""
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
M26 = HERE.parent
MSDB = Path.home() / "Documents/ModelSEEDDatabase"
KEGG_FILE = HERE / "kegg_cache/link_reaction_ko_2026-11-20.tsv"
OUT = HERE / "maps/term_maps.json"
TEMPLATES = ["GramNegModelTemplateV7.json", "GramPosModelTemplateV7.json", "ArchaeaTemplateV6.json", "Core-V6.json"]
EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def git(*a):
    return subprocess.run(["git", "-C", str(MSDB), *a], capture_output=True, text=True, check=True).stdout.strip()


def main():
    # KEGG KO -> R
    ko2r = defaultdict(set)
    for line in KEGG_FILE.read_text().splitlines():
        ko, rn = line.split("\t")
        ko2r[ko.replace("ko:", "")].add(rn.replace("rn:", ""))
    # ModelSEED aliases
    al = MSDB / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt"
    ecf = MSDB / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt"
    r2ms, bigg2ms = defaultdict(set), defaultdict(set)
    for line in al.read_text().splitlines()[1:]:
        rid, ext, src = line.split("\t")
        if src == "KEGG":
            r2ms[ext].add(rid)
        elif src == "BiGG":
            bigg2ms[ext].add(rid)
    ec2ms = defaultdict(set)
    for line in ecf.read_text().splitlines()[1:]:
        rid, ec, src = line.split("\t")
        if EC_RE.match(ec):
            ec2ms[ec].add(rid)
    ko2ms = {}
    for ko, rs in ko2r.items():
        s = set().union(*(r2ms.get(r, set()) for r in rs))
        if s:
            ko2ms[ko] = sorted(s)
    # coverage in templates
    cov = {}
    for t in TEMPLATES:
        tj = json.load(open(M26 / "templates" / t))
        trx = {r["id"].rsplit("_", 1)[0] for r in tj["reactions"]}
        cov[t] = {
            "template_reactions": len(trx),
            "ko_with_template_rxn": sum(1 for v in ko2ms.values() if trx & set(v)),
            "ec_with_template_rxn": sum(1 for v in ec2ms.values() if trx & set(v)),
            "bigg_with_template_rxn": sum(1 for v in bigg2ms.values() if trx & v),
        }
    out = {
        "provenance": {
            "kegg_link_file": str(KEGG_FILE.relative_to(M26)), "kegg_link_md5": md5(KEGG_FILE),
            "kegg_link_source": "https://rest.kegg.jp/link/reaction/ko", "kegg_link_downloaded": "2026-11-20T21:09Z",
            "kegg_link_rows": sum(len(v) for v in ko2r.values()), "kegg_kos": len(ko2r),
            "modelseed_db_path": str(MSDB), "modelseed_db_branch": git("branch", "--show-current"),
            "modelseed_db_commit": git("rev-parse", "HEAD"),
            "reaction_aliases_md5": md5(al), "reaction_ecs_md5": md5(ecf),
            "kegg_R_with_modelseed_alias": len(r2ms), "ec_rule": "complete 4-level EC numbers only",
        },
        "counts": {"ko2rxn": len(ko2ms), "ec2rxn": len(ec2ms), "bigg2rxn": len(bigg2ms)},
        "template_coverage": cov,
        "ko2rxn": ko2ms,
        "ec2rxn": {k: sorted(v) for k, v in ec2ms.items()},
        "bigg2rxn": {k: sorted(v) for k, v in bigg2ms.items()},
    }
    OUT.parent.mkdir(exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=1, sort_keys=False)
    print(json.dumps({k: out[k] for k in ("provenance", "counts", "template_coverage")}, indent=1))


if __name__ == "__main__":
    main()
