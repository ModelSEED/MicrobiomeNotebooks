"""Shared paths/helpers for the biology/annotation check. READ-ONLY on the salterns repo."""
import csv, json, re, sys
from pathlib import Path
sys.dont_write_bytecode = True
S = Path("/home/freiburger/Documents/MicrobiomeNotebooks/salterns")
VB = Path("/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/a24fbb2e-9d33-5e05-aeea-a07d777ccde7/scratchpad/mepn_verify/verify_bio")
MEPN_MODELS = S / "models" / "MePn_consumers"
RAST_DIR = S / "jplfaria" / "Cliff MAGs" / "RASTified_MAG_json"
GUILDS = S / "correlational" / "r232" / "guilds_r232.csv"
IDS = json.load(open(S / "correlational" / "iterativeIDs_r232.json"))
OBLIGATE = {"Roseovarius.9", "Roseovarius.11", "Albimonas.1"}

def model_ids():
    return sorted(p.name[:-5] for p in MEPN_MODELS.glob("*.json"))

def mag_of(mid):
    return mid.replace(".contigs", "")

def iid(mid):
    return IDS[mag_of(mid)]

def guilds():
    with open(GUILDS) as fh:
        return {r["mag"]: r for r in csv.DictReader(fh)}

def rast(mid):
    return json.load(open(RAST_DIR / f"{mid}__.RAST.json"))

def short(fid):
    return fid.split("RAST.")[-1]
