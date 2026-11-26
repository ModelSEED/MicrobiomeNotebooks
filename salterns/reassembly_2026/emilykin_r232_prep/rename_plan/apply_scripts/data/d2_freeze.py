"""D2: freeze the R220-era network MAG link BEFORE the R232 map is written.
Writes (in ~/Documents/EmilyKin):
  network/mag_linked_asvs_r220_frozen.tsv          asv_iterativeID, MAG  (140 rows; LF)
  network/mag_guild_iterative_ID_stale_frozen.tsv   MAG, stale_iterative_ID (291 rows; LF) -- the guild table's
      iterative_ID column at 9e087c3 (6a04588-era map), read only by the root spieceasi_modules_flanking.py audit.
Checks: every frozen ASV is a network node (abundances.csv column); the ASV's asv_to_mag_mapping best_mag == bin;
every consumer universe (SPIEC-EASI adjacency/edge nodes, module members) is a subset of abundances.csv columns,
so frozen-set membership == HEAD-map membership for every node any consumer can see."""
import csv, json, glob, os, sys, subprocess
R = os.path.expanduser("~/Documents/EmilyKin")
head = subprocess.run(["git", "-C", R, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
assert head == "9e087c3", head
assert subprocess.run(["git", "-C", R, "status", "--porcelain", "--", "mag_iterativeID_old_to_new.json",
                       "mag_guild_classification.csv"], capture_output=True, text=True).stdout == "", "map/guild already modified"
bin2id = json.load(open(f"{R}/mag_iterativeID_old_to_new.json"))
id2bin = {v: k for k, v in bin2id.items()}
assert len(id2bin) == len(bin2id) == 291
cols = open(f"{R}/abundances.csv").readline().rstrip("\r\n").split(",")[1:]
assert len(cols) == 935, len(cols)
frozen = [(a, id2bin[a]) for a in cols if a in id2bin]
print("frozen pairs:", len(frozen))
assert len(frozen) == 140
a2m = json.load(open(f"{R}/asv_to_mag_mapping.json"))["asv_to_mag"]
ids = json.load(open(f"{R}/iterativeIDs.json"))           # ASV iterativeID -> hash (as in recompute_mag_ids.py)
id2hash = dict(ids)
bad = [(a, b, a2m.get(id2hash.get(a), {}).get("best_mag")) for a, b in frozen if a2m.get(id2hash.get(a), {}).get("best_mag") != b]
print("best_mag != bin:", bad)
assert not bad
# consumer universes
U = set(cols)
uni = {}
for f in glob.glob(f"{R}/network/network_module_membership_p_value_FDR*.json"):
    d = json.load(open(f))
    uni[os.path.relpath(f, R)] = {m for k, v in d.items() if isinstance(v, dict) for m in v.get("members", [])}
for f in glob.glob(f"{R}/network/compositional_networks/out/asv/*.csv"):
    rows = list(csv.reader(open(f)))
    hdr = rows[0]
    s = set()
    if "source" in hdr and "target" in hdr:
        i, j = hdr.index("source"), hdr.index("target"); s = {r[i] for r in rows[1:]} | {r[j] for r in rows[1:]}
    elif hdr and hdr[0] == "":     # square matrix (adjacency / weights)
        s = set(hdr[1:])
    elif "node" in hdr:
        s = {r[hdr.index("node")] for r in rows[1:]}
    if s:
        uni[os.path.relpath(f, R)] = s
for k, s in sorted(uni.items()):
    extra = s - U
    print(f"universe {k}: {len(s)} nodes, not in abundances.csv: {len(extra)}")
    assert not extra, (k, sorted(extra)[:5])
out = f"{R}/network/mag_linked_asvs_r220_frozen.tsv"
assert not os.path.exists(out)
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["asv_iterativeID", "MAG"])
    for a, b in sorted(frozen):
        w.writerow([a, b])
# stale guild-table column (root spieceasi_modules_flanking.py audit only)
g = list(csv.DictReader(open(f"{R}/mag_guild_classification.csv")))
assert len(g) == 291
out2 = f"{R}/network/mag_guild_iterative_ID_stale_frozen.tsv"
assert not os.path.exists(out2)
with open(out2, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["MAG", "stale_iterative_ID"])
    for r in g:
        w.writerow([r["MAG"], r["iterative_ID"]])
print("wrote", out, "and", out2)
