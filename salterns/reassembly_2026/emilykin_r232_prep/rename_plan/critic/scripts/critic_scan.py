#!/usr/bin/env python
"""Independent completeness scan for the EmilyKin R220 -> R232 MAG rename (critic part).

Read-only on ~/Documents/EmilyKin. Scans tracked files AND untracked/ignored files (outside the
raw-read / synthetic-genome trees) for the 51 old MAG IDs and the 4 ASV-shared new IDs, in
strict form and in variant encodings that a boundary-exact tokenizer misses:
  strict      ID exactly, boundaries (?<![A-Za-z0-9])ID(?![0-9A-Za-z]) (leading '_' allowed -> 'embedded')
  ca_display  'Ca. X.N', 'Candidatus X.N', bare 'X.N' for Ca_X roots
  sep         root + '_' / '-' / ' ' + number (e.g. Ignavibacteria_5, Giesbergeria_1_m, file-safe names)
  escaped     markdown/LaTeX escaped underscores (Azospira\\_A.1\\_m)
  split_run   docx/pptx: ID only visible after joining the runs of a paragraph
Container formats: .gz, xlsx/docx/pptx (all XML parts incl. charts, comments, embedded OOXML),
pdf (pdftotext), FASTA (header lines only). Everything else is scanned as bytes.
Output: critic/scan/hits.tsv (one row per file x id x variant, count + first context).
"""
import gzip, io, json, os, re, subprocess, sys, zipfile, csv
from collections import defaultdict

REPO = os.path.expanduser('~/Documents/EmilyKin')
PLAN = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # rename_plan/
OUT = os.path.join(PLAN, 'critic', 'scan')
os.makedirs(OUT, exist_ok=True)
CAP = 400 * 1024 * 1024

ch = list(csv.DictReader(open(os.path.join(PLAN, 'data', 'mag_id_changes_r220_r232.csv'))))
OLD = {r['old_mag_iterativeID']: r['MAG'] for r in ch}
NEW = {r['new_mag_iterativeID']: r['MAG'] for r in ch}
NEW_ASV = [r['new_mag_iterativeID'] for r in ch if r['new_is_asv_id'] == 'True']
IDS = {}
for i in OLD: IDS[i] = 'old'
for i in NEW_ASV: IDS[i] = 'new_asv_shared'


def split_id(i):
    m = re.match(r'^(.*)\.(\d+)(_m)?$', i)
    return m.group(1), m.group(2), bool(m.group(3))


def esc_root(root, flex_us):
    parts = [re.escape(p) for p in root.split('_')]
    us = r'(?:_|\\_|\\\\_)' if flex_us else '_'
    return us.join(parts)


PATS = []  # (variant, compiled bytes regex with named groups g<k> -> IDLIST[k])
IDLIST = list(IDS)
_alts = defaultdict(list)
for k, i in enumerate(IDLIST):
    root, n, minted = split_id(i)
    tail_m = r'_m(?![A-Za-z0-9])' if minted else r'(?![0-9])(?!(?:\\?_)?m(?![a-z]))(?![A-Za-z])'
    tail_m_flex = r'(?:_|\\_|\\\\_)m(?![A-Za-z0-9])' if minted else tail_m
    strict_tail = r'(?![0-9A-Za-z])' if minted else r'(?![0-9])(?!_m)(?![A-Za-z])'
    _alts['strict'].append(f'(?P<g{k}>' + re.escape(i) + ')' + strict_tail)
    base_root = root[3:] if root.startswith('Ca_') else None
    if base_root:
        _alts['ca_display'].append(f'(?P<g{k}>(?<![A-Za-z0-9_])(?:Ca\\.\\s?|Ca\\s|Candidatus\\s|\\*Ca\\.\\*\\s?|)' + re.escape(base_root) + r'\.' + n + ')' + tail_m)
    _alts['sep'].append(f'(?P<g{k}>' + esc_root(root, False) + r'[_\- ]' + n + ('_m' if minted else '') + ')' + strict_tail)
    _alts['escaped_any'].append(f'(?P<g{k}>' + esc_root(root, True) + r'\.' + n + (r'(?:_|\\_|\\\\_)m' if minted else '') + ')' + strict_tail)
for var, alts in _alts.items():
    PATS.append((var, re.compile((r'(?<![A-Za-z0-9])(?:' + '|'.join(alts) + ')').encode())))

ROOTS = sorted({split_id(i)[0].replace('Ca_', '') for i in IDS} | {split_id(i)[0].split('_')[0] for i in IDS})
ROOTS = [r for r in ROOTS if r != 'Ca']
ROOT_BYTES = [r.encode() for r in ROOTS]

FASTA_EXT = {'.fa', '.faa', '.fna', '.fasta', '.ffn', '.fas', '.aln', '.afa', '.ffa'}
SKIP_DIRS = ('meta/', 'mag/', 'synthetic_genomes', 'kofamscan/', 'bakta/', 'MAG_faas/', '.claude/', 'kofam_new/profiles')
SKIP_EXT = {'.mmi', '.bam', '.sam', '.fastq', '.png', '.jpg', '.jpeg', '.gif', '.tif', '.tiff', '.hmm', '.zip', '.tar'}


def files():
    tr = subprocess.run(['git', '-C', REPO, 'ls-files', '-z'], capture_output=True).stdout.split(b'\0')
    un = subprocess.run(['git', '-C', REPO, 'ls-files', '-o', '-z'], capture_output=True).stdout.split(b'\0')
    out = []
    for status, lst in (('tracked', tr), ('untracked', un)):
        for p in lst:
            if not p:
                continue
            p = p.decode('utf-8', 'surrogateescape')
            if p.startswith(SKIP_DIRS) or '/._' in p or p.startswith('._') or '__pycache__' in p or p.startswith('bins_'):
                continue
            if any(p.endswith(s) for s in ('.fastq.gz', '.fq.gz', '.tar.gz')):
                continue
            if status == 'untracked' and p.startswith(UNTRACKED_SKIP):
                SKIPPED.append((status, p, 'untracked ASV-run tree'))
                continue
            if not os.path.isfile(os.path.join(REPO, p)):
                SKIPPED.append((status, p, 'not a regular file (directory / nested repo)'))
                continue
            try:
                sz = os.path.getsize(os.path.join(REPO, p))
            except OSError:
                continue
            ext = os.path.splitext(p)[1].lower()
            if sz > BIG and ext not in FASTA_EXT:
                SKIPPED.append((status, p, f'{sz // 1000000} MB > cap'))
                continue
            out.append((status, p))
    return out


UNTRACKED_SKIP = ('community_modeling/runs/', 'community_modeling/_pinned/', 'community_modeling/plan_evidence/',
                  'metabolic_complementarity/profiles/', 'community_modeling/analysis/data_inconsistencies/')
BIG = 60 * 1000 * 1000
SKIPPED = []


def ooxml_texts(data, name):
    """yield (part, text) for every XML part; for docx/pptx also paragraph-joined text."""
    try:
        z = zipfile.ZipFile(io.BytesIO(data))
    except Exception:
        return
    for n in z.namelist():
        b = z.read(n)
        if n.lower().endswith(('.xlsx', '.docx', '.pptx', '.xlsm')):
            for part, t in ooxml_texts(b, n):
                yield n + '>' + part, t
            continue
        if not n.endswith(('.xml', '.rels', '.vml')):
            continue
        yield n, b
        if re.search(rb'<w:p[ >]|<a:p>', b):
            paras = re.split(rb'</w:p>|</a:p>', b)
            joined = []
            for para in paras:
                ts = re.findall(rb'<(?:w|a):t(?: [^>]*)?>([^<]*)</(?:w|a):t>', para)
                if ts:
                    joined.append(b''.join(ts))
            yield n + '#joined', b'\n'.join(joined)


def blobs(path):
    full = os.path.join(REPO, path)
    try:
        sz = os.path.getsize(full)
    except OSError:
        return
    ext = os.path.splitext(path)[1].lower()
    if ext in SKIP_EXT or sz > CAP:
        if ext in FASTA_EXT or path.endswith('.fa.gz'):
            pass
        else:
            return
    if ext in FASTA_EXT:
        r = subprocess.run(['grep', '-a', '^>', full], capture_output=True)
        yield 'headers', r.stdout
        return
    data = open(full, 'rb').read()
    if ext == '.gz':
        try:
            data = gzip.decompress(data)
        except Exception:
            pass
        yield 'gunzip', data
        return
    if ext in ('.xlsx', '.docx', '.pptx', '.xlsm', '.qza', '.qzv', '.npz'):
        if ext in ('.qza', '.qzv', '.npz'):
            try:
                z = zipfile.ZipFile(io.BytesIO(data))
                for n in z.namelist():
                    yield n, z.read(n)
            except Exception:
                yield 'raw', data
            return
        for part, t in ooxml_texts(data, path):
            yield part, t
        return
    if ext == '.pdf':
        r = subprocess.run(['pdftotext', '-layout', full, '-'], capture_output=True)
        yield 'pdftotext', r.stdout
        return
    yield 'raw', data


def main():
    fl = files()
    print('files', len(fl), file=sys.stderr)
    rows = []
    pathhits = []
    for k, (status, p) in enumerate(fl):
        for i in IDS:
            if re.search(r'(?<![A-Za-z0-9])' + re.escape(i) + r'(?![0-9A-Za-z])', p):
                pathhits.append((status, p, i))
        counts = defaultdict(int)
        ctx = {}
        for part, b in blobs(p):
            if not b:
                continue
            # fast literal search for roots -> merged windows; precise regexes run only in windows
            wins = []
            for r in ROOT_BYTES:
                j = b.find(r)
                while j != -1:
                    wins.append((max(0, j - 40), j + len(r) + 30))
                    j = b.find(r, j + 1)
            if not wins:
                continue
            wins.sort()
            merged = [list(wins[0])]
            for a, e in wins[1:]:
                if a <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([a, e])
            joined = part.endswith('#joined')
            strict_spans = defaultdict(set)
            seen = set()
            for var, rx in PATS:
                for a, e in merged:
                    for m in rx.finditer(b, a, min(len(b), e)):
                        if (var, m.start()) in seen:
                            continue
                        seen.add((var, m.start()))
                        i = IDLIST[int(m.lastgroup[1:])]
                        v = var
                        if var == 'strict':
                            pre = b[max(0, m.start() - 1):m.start()]
                            if pre == b'_':
                                v = 'strict_embedded'
                            strict_spans[i].add(m.start())
                        else:
                            if any(abs(m.start() - s0) <= 12 for s0 in strict_spans[i]) or m.group(0) == i.encode():
                                continue
                        key = (part.split('#')[0] if not joined else 'joined', i, v)
                        counts[key] += 1
                        if key not in ctx:
                            ctx[key] = b[max(0, m.start() - 80):m.end() + 60].decode('utf-8', 'replace').replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')
        # collapse: for docx/pptx the joined text duplicates raw xml hits; keep raw xml hits, and
        # mark joined-only ids as split_run
        per_id_raw = defaultdict(int)
        for (part, i, v), c in counts.items():
            if part != 'joined' and v.startswith('strict'):
                per_id_raw[i] += c
        for (part, i, v), c in counts.items():
            if part == 'joined':
                if v.startswith('strict'):
                    extra = c - per_id_raw.get(i, 0)
                    if extra > 0:
                        rows.append((status, p, i, IDS[i], 'split_run', extra, ctx[(part, i, v)]))
                continue
            rows.append((status, p, i, IDS[i], v, c, ctx[(part, i, v)] + ('' if part in ('raw',) else f'  [part={part}]')))
        if k % 2000 == 0:
            print(k, len(rows), file=sys.stderr)
    with open(os.path.join(OUT, 'hits.tsv'), 'w') as f:
        f.write('status\tfile\tid\tid_kind\tvariant\tcount\tcontext\n')
        for r in rows:
            f.write('\t'.join(map(str, r)) + '\n')
    with open(os.path.join(OUT, 'path_hits.tsv'), 'w') as f:
        f.write('status\tfile\tid\n')
        for r in pathhits:
            f.write('\t'.join(r) + '\n')
    with open(os.path.join(OUT, 'skipped.tsv'), 'w') as f:
        f.write('status\tfile\treason\n')
        for r in SKIPPED:
            f.write('\t'.join(r) + '\n')
    print('rows', len(rows), 'pathhits', len(pathhits), 'skipped', len(SKIPPED), file=sys.stderr)


if __name__ == '__main__':
    main()
