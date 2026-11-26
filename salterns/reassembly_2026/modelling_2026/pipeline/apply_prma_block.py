#!/usr/bin/env python
"""Close rxn00844 (annotated by RAST as 'methane monooxygenase component A alpha chain', EC 1.14.13.25) in the two
Roseovarius MePn-consumer models that carry it. The gene is a propane 2-monooxygenase (PrmA): KOfam PrmA score 1016
(threshold 744) vs MmoX 228 (threshold 1169); 65 % identity to Gordonia/Rhodococcus PrmA, 32-34 % to MmoX; prmABCD operon
with GroEL and an alcohol dehydrogenase; no mmo/pmo genes in the genomes (models/methanogen/CURATION.md, Roseovarius.2
section). Without the block the community protocol oxidises the MePn-derived methane through it.
Run once: ~/Documents/py_venv/bin/python apply_prma_block.py   (idempotent; rewrites models/MANIFEST.md5 entries)
"""
import hashlib
import json
from pathlib import Path

import cobra

M = Path(__file__).resolve().parents[1] / "models"
TARGETS = ["R2_restored_C_black__metabat2.bin.011_sub", "R2_B_D2__metabat2.bin.006_sub"]
for mag in TARGETS:
    f = M / "mepn_consumer" / f"{mag}.json"
    m = cobra.io.load_json_model(str(f))
    r = m.reactions.get_by_id("rxn00844_c0")
    r.bounds = (0.0, 0.0)
    r.notes["salterns_2026_block"] = "closed 2026-11-24: gene is propane 2-monooxygenase (PrmA), not methane monooxygenase"
    cobra.io.save_json_model(m, str(f))
    print(mag, "rxn00844_c0 closed;", hashlib.md5(f.read_bytes()).hexdigest())
