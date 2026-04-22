"""Beautify methylphosphonate_map.json and betaine_map.json using editEscher.

Pipeline:
1. ``cleanEscherJSON`` — abbreviate metabolite + reaction labels. Metabolite
   abbreviations come from ``compound_abbrev_mapping.json`` (BiGG whenever
   available; proposed for the four methylphosphonate intermediates that
   lack BiGG/MetaCyc/AraCyc short forms).
2. ``apply_layout`` — apply the manual label positions fine-tuned in the
   *_edited.json* files, promote the four highlighted compounds (betaine,
   methylphosphonate, TMA, methane) to primary nodes so Escher draws them
   larger, and drop leftover custom text labels (MePn map only).
3. ``reverse_reactions_in_map`` (betaine only) — flip rxn17220 so its
   arrows point down.
"""
import json
import sys
from pathlib import Path
sys.path.insert(0, '/Users/andrewfreiburger/Documents/Research/editEscher/src')

from escher_edit import (
    cleanEscherJSON,
    apply_layout,
    reverse_reactions_in_map,
    EscherSVG_processing,
    EscherStyle,
)


# ---- Metabolite abbreviations (BiGG-first, sourced from
# compound_abbrev_mapping.json — see that file for {abbrev, source}).
with open('compound_abbrev_mapping.json') as _fh:
    _COMPOUND_MAPPING = json.load(_fh)

# The methylphosphonate map uses compartment-suffixed bigg_ids (``cpd00001_c0``),
# while the betaine map uses bare IDs (``cpd00001``). Build both variants.
_MAP_BY_CPD = {cid: v['abbrev'] for cid, v in _COMPOUND_MAPPING.items()}

MEPN_CPD_ABBREV = {f'{cid}_c0': abbr for cid, abbr in _MAP_BY_CPD.items()}
MEPN_CPD_ABBREV['cpd01024_e0'] = _MAP_BY_CPD['cpd01024']  # methane is extracellular in MePn map

BET_CPD_ABBREV = dict(_MAP_BY_CPD)


MEPN_RXN_ABBREV = {
    'rxs43120_c0': 'C-P lyase',
    'rxn26457_c0': 'RPnTP hydrolase',
    'rxn26456_c0': 'PRPn transferase',
}

BET_RXN_ABBREV = {
    'rxn17220': 'Betaine reductase',
    'rxn10613': 'TMA methyltrfrs',
    'rxn10496': 'DMA methyltrfrs',
    'rxn10572': 'MMA methyltrfrs',
    'rxn10569': 'McrA',
}


# ---- Manual label positions (extracted verbatim from *_edited.json) ----

MEPN_NODE_LABEL_POSITIONS = {
    '22': (-32.26, 752.08),    # NADH
    '26': (188.23, 1014.71),   # NAD
    '27': (146.84, 769.97),    # amet (SAM)
    '28': (260.27, 976.69),    # met__L (Methionine)
    '29': (117.31, 1056.23),   # ch4
    '30': (-49.02, 1005.81),   # dad_5 (dAdo)
    '31': (-124.74, 958.11),   # prpcp
    '33': (189.28, 405.01),    # h2o
    '37': (187.32, 615.04),    # ppi
    '38': (-11.69, 600.18),    # h
    '40': (118.18, 720.31),    # prpn
    '41': (206.14, 64.61),     # atp
    '45': (207.23, 241.23),    # ade
    '46': (128.28, 14.95),     # mepn
    '47': (124.18, 334.65),    # rpntp
}

MEPN_RXN_LABEL_POSITIONS_BY_ID = {
    '3': (134.30, 894.80),   # C-P lyase
    '4': (122.30, 512.90),   # RPnTP hydrolase
    '5': (116.30, 155.10),   # PRPn transferase
}

BET_NODE_LABEL_POSITIONS = {
    '0':  (-263.71,  -24.68),  # trdox
    '4':  (-259.70, -233.01),  # trdrd
    '5':  (-139.40, -272.32),  # glyb (BET)
    '6':  (  -8.19, -240.04),  # h
    '7':  ( -85.13,   46.37),  # tma
    '8':  ( -20.74,  -26.27),  # h2o
    '9':  (  -9.88,  -70.65),  # actp
    '10': (   9.43, -191.82),  # pi
    '14': (-252.55,   83.27),  # com
    '15': ( -79.11,  299.14),  # dma
    '16': (-268.53,  199.70),  # mcom
    '20': (-256.11,  337.51),  # com
    '21': ( -86.21,  567.99),  # mma
    '22': (-269.42,  473.88),  # mcom
    '26': (-248.57,  594.36),  # com
    '27': ( -14.93,  745.61),  # nh3
    '28': (-254.78,  785.73),  # mcom
    '32': ( -21.16,  992.18),  # hsfd
    '33': (-152.54, 1072.10),  # ch4
    '34': (-236.61,  993.29),  # h
    '35': (-236.31,  825.57),  # cob (HTP)
}

BET_RXN_LABEL_POSITIONS_BY_ID = {
    '0': ( -89.92, -120.23),  # Betaine reductase
    '1': ( -97.98,  172.31),  # TMA methyltrfrs
    '2': ( -97.98,  430.49),  # DMA methyltrfrs
    '3': ( -97.69,  679.07),  # MMA methyltrfrs
    '4': (-100.36,  920.33),  # McrA
}


# Nodes to highlight (promote to primary so they render larger).
# After cleanEscherJSON the bigg_ids match the new abbreviations.
MEPN_HIGHLIGHT = (_MAP_BY_CPD['cpd25960'], _MAP_BY_CPD['cpd01024'])  # mepn, ch4
BET_HIGHLIGHT = (_MAP_BY_CPD['cpd00540'], _MAP_BY_CPD['cpd00441'], _MAP_BY_CPD['cpd01024'])  # glyb, tma, ch4


if __name__ == "__main__":
    # --- methylphosphonate map ---
    cleanEscherJSON(
        "methylphosphonate_map.json",
        abbrev_map=MEPN_CPD_ABBREV,
        rxn_abbrev_map=MEPN_RXN_ABBREV,
    )
    apply_layout(
        "methylphosphonate_map_cleaned0.json",
        output_path="methylphosphonate_map_cleaned.json",
        metabolite_label_dx=None,
        metabolite_label_dy=None,
        node_label_positions=MEPN_NODE_LABEL_POSITIONS,
        reaction_label_positions_by_id=MEPN_RXN_LABEL_POSITIONS_BY_ID,
        primary_bigg_ids=MEPN_HIGHLIGHT,
        demote_other_primaries=True,
        highlight_label_size=30,
        highlight_offset_scale=1.8,
        remove_all_text_labels=True,
    )

    # --- betaine map ---
    cleanEscherJSON(
        "betaine_map.json",
        abbrev_map=BET_CPD_ABBREV,
        rxn_abbrev_map=BET_RXN_ABBREV,
    )
    apply_layout(
        "betaine_map_cleaned0.json",
        output_path="betaine_map_layout.json",
        metabolite_label_dx=None,
        metabolite_label_dy=None,
        node_label_positions=BET_NODE_LABEL_POSITIONS,
        reaction_label_positions_by_id=BET_RXN_LABEL_POSITIONS_BY_ID,
        primary_bigg_ids=BET_HIGHLIGHT,
        demote_other_primaries=True,
        highlight_label_size=30,
        highlight_offset_scale=1.8,
    )
    reverse_reactions_in_map(
        "betaine_map_layout.json",
        rxn_bigg_ids=("Betaine reductase",),
        output_path="betaine_map_cleaned.json",
    )

    # --- Optional SVG post-processing ---
    # Escher's Python Builder only outputs HTML. To get an SVG, export from
    # the Escher web UI (File > Export as SVG) after loading the cleaned
    # JSON, and save it as ``<stem>_cleaned.svg`` alongside the JSON. The
    # block below picks up whatever SVGs exist and enlarges the highlighted
    # node labels by reading the ``label_size`` field from the JSON.
    svg_style = EscherStyle(node_label_px=10, rxn_label_px=30)
    for json_path, svg_path in [
        ("methylphosphonate_map_cleaned.json", "methylphosphonate_map_cleaned.svg"),
        ("betaine_map_cleaned.json", "betaine_map_cleaned.svg"),
    ]:
        if Path(svg_path).exists():
            EscherSVG_processing(
                svg_path,
                style=svg_style,
                json_path=json_path,  # reads per-node label_size values
            )
            print(f"Post-processed: {svg_path}")
        else:
            print(f"Skipped (no SVG): {svg_path} — export it from Escher UI to enable font enlargement")
