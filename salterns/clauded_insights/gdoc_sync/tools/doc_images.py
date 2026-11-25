#!/usr/bin/env python
"""Images for the Google Doc figure swaps (figure_swap.py).

Google's insertInlineImage fetches the image from a public URL and refuses images above 25 megapixels, so figures larger
than MAX_MP are downscaled (Lanczos) into gdoc_sync/doc_images/; smaller ones are used as they are. The images are
fetched from raw.githubusercontent.com at a commit that holds them, so commit and push before running figure_swap.py.
Writes gdoc_sync/doc_images/images.json: figure -> {path (relative to salterns/), width_px, height_px}.
Run: ~/Documents/py_venv/bin/python doc_images.py
"""
import json
from pathlib import Path

from PIL import Image

Image.MAX_IMAGE_PIXELS = None
S = Path(__file__).resolve().parents[3]                 # salterns/
OUT = S / "clauded_insights" / "gdoc_sync" / "doc_images"
MAX_MP = 20.0
FIGURES = {
    "Figure 1": "reassembly_2026/community_ko/figures/figure1_2026.png",
    "Figure 7": "reassembly_2026/modelling_2026/simulations/figure7_pergram.png",
    "Figure 8": "correlational/final585/within_condition_louvain_network_final585.png",
    "Figure 9": "correlational/final585/guild_methane_raw_final585.png",
    "Figure 10": "reassembly_2026/modelling_2026/synthesis/figure10_pathways.png",
    "Figure S4": "reassembly_2026/community_ko/figures/figureS4_2026.png",
    "Figure S6": "reassembly_2026/community_ko/figures/figureS6_2026.png",
    "Figure S8": "correlational/final585/guild_methane_dnanorm_final585.png",
    "Figure S10": "reassembly_2026/modelling_2026/figureS10/figureS10_final.png",
    "Figure S11": "correlational/final585/figureS11_final585.png",
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out = {}
    for fig, rel in FIGURES.items():
        src = S / rel
        if not src.exists():
            print(f"{fig}: {rel} not built yet, skipped")
            continue
        im = Image.open(src)
        w, h = im.size
        if w * h / 1e6 > MAX_MP:
            f = (MAX_MP * 1e6 / (w * h)) ** 0.5
            im = im.convert("RGBA").resize((int(w * f), int(h * f)), Image.LANCZOS)
            dst = OUT / (fig.replace(" ", "") + "_" + src.name)
            im.save(dst, optimize=True)
            rel = str(dst.relative_to(S))
            w, h = im.size
        out[fig] = {"path": rel, "width_px": w, "height_px": h}
        print(f"{fig}: {rel} {w}x{h}")
    json.dump(out, open(OUT / "images.json", "w"), indent=1)


if __name__ == "__main__":
    main()
