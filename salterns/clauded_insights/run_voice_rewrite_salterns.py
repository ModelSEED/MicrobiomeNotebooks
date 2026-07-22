"""
Restyle the CORRECTED saltern manuscript (manuscript_corrected.md) into Andrew
Freiburger's writing voice with the local fine-tuned model (`andrew-voice`,
Mistral-NeMo LoRA served by ollama), PARAGRAPH BY PARAGRAPH so the model cannot
summarise across paragraphs.

Safeguards (a manuscript must not lose facts):
  * headings, figure/table captions (**...**), and the References / Figures and Tables
    / Data Availability / Acknowledgments sections are passed through verbatim;
  * every prose paragraph is rewritten alone with a strict "do not shorten / preserve
    all citations & numbers" prompt;
  * a CITATION GUARD compares the numeric-citation set before/after -- if the rewrite
    drops (or invents) any citation number, the ORIGINAL paragraph is kept.
Output: manuscript_corrected.andrew_voice.md  (a stylistic draft -- still fact-check).
stdlib only.  ~/Documents/py_venv/bin/python.
"""
from __future__ import annotations
import json, re, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "manuscript_corrected.md"
OUT = HERE / "manuscript_corrected.andrew_voice.md"
API = "http://localhost:11434/api/generate"
MODEL = "andrew-voice"
PASSTHROUGH = re.compile(r"^# (References|Figures and Tables|Data Availability|Acknowledgments)")
CITE = re.compile(r"\((\d+(?:\s*[,–-]\s*\d+)*)\)")   # (34), (41–44), (2, 5)

PROMPT = (
    "Rewrite this SINGLE paragraph from a research paper in Andrew Freiburger's academic voice "
    "(precise, scholarly, fluent). Rules: do NOT summarise, shorten, merge, or omit any content; "
    "keep it about the same length; preserve EVERY parenthetical citation such as (34) or (41-44), "
    "every number, statistic (R2, rho, p, %), gene, organism, and figure/table reference EXACTLY; "
    "do not invent citations or facts. Output ONLY the rewritten paragraph.\n\n"
    "PARAGRAPH:\n{chunk}\n\nREWRITTEN PARAGRAPH:\n"
)


def cites(s: str) -> set[str]:
    out = set()
    for grp in CITE.findall(s):
        for tok in re.split(r"[,\-–]", grp):
            tok = tok.strip()
            if tok.isdigit():
                out.add(tok)
    return out


def gen(prompt: str, npred: int) -> str:
    body = json.dumps({"model": MODEL, "prompt": prompt, "stream": False, "keep_alive": "15m",
                       "options": {"num_ctx": 8192, "num_predict": npred, "temperature": 0.55,
                                   "top_p": 0.9, "repeat_penalty": 1.1}}).encode()
    req = urllib.request.Request(API, data=body, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=900) as r:
        return json.loads(r.read())["response"]


def main():
    text = SRC.read_text()
    m = re.search(r"(?m)^# Abstract", text)
    header, body = text[:m.start()].rstrip(), text[m.start():]
    paras = [p for p in body.split("\n") if p.strip()]   # one line == one heading/caption/paragraph

    out, passthrough, kept, rew, guarded = [header], False, 0, 0, 0
    for i, p in enumerate(paras, 1):
        first = p.lstrip().splitlines()[0]
        if first.startswith("# "):
            passthrough = bool(PASSTHROUGH.match(first))
        # verbatim: any heading line, caption paragraphs, or inside a passthrough section
        if passthrough or first.startswith("#") or p.lstrip().startswith("**"):
            out.append(p); kept += 1
            OUT.write_text("\n\n".join(out) + "\n"); continue
        words = len(p.split())
        try:
            resp = gen(PROMPT.format(chunk=p), min(3072, max(256, int(words * 2.6))))
            resp = re.split(r"\n(?:REWRITTEN PARAGRAPH:|PARAGRAPH:)", resp)[0].strip()
        except Exception as e:
            resp = ""; print(f"    !! {e}", flush=True)
        # citation guard: keep original if any citation dropped or invented, or empty/too short
        if not resp or cites(resp) != cites(p) or len(resp.split()) < 0.6 * words:
            out.append(p); guarded += 1
            print(f"[{i}/{len(paras)}] GUARD-kept original ({words}w; cites {sorted(cites(p))})", flush=True)
        else:
            out.append(resp); rew += 1
            print(f"[{i}/{len(paras)}] rewrote ({words}w -> {len(resp.split())}w)", flush=True)
        OUT.write_text("\n\n".join(out) + "\n")

    OUT.write_text("\n\n".join(out) + "\n")
    print(f"\nDONE  rewrote={rew}  guard-kept={guarded}  verbatim={kept}  -> {OUT.name}", flush=True)


if __name__ == "__main__":
    main()
