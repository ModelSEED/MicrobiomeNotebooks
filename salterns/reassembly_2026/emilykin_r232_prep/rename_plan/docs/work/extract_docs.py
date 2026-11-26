"""Read-only text extraction of EmilyKin documents into segments.jsonl.
docx: every w:p in document/footnotes/endnotes/comments/headers/footers, with final (accepted) text and deleted text
separately; tables keep (table,row,cell) + whole-row text. xlsx: every non-empty cell with sheet, header, row bins.
pptx: shape/table/notes text per slide. md/html: lines. pdf: pdftotext -layout pages."""
import json, re, subprocess, sys, zipfile
from pathlib import Path
from lxml import etree
EK = Path.home()/'Documents/EmilyKin'
OUT = Path(__file__).parent/'segments.jsonl'
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
BIN = re.compile(r'\b(?:CAN_\d_bin\.\d+|coasm_bin\.\d+)\b')

def ptext(p):
    fin, dele = [], []
    for el in p.iter(W+'t', W+'delText', W+'tab', W+'br'):
        if el.tag == W+'t':
            fin.append(el.text or '')
        elif el.tag == W+'delText':
            dele.append(el.text or '')
        else:
            fin.append(' ')
    return ''.join(fin), ''.join(dele)

def docx_segments(path):
    z = zipfile.ZipFile(path)
    parts = [n for n in z.namelist() if re.match(r'word/(document|footnotes|endnotes|comments|header\d*|footer\d*)\.xml$', n)]
    for part in parts:
        root = etree.fromstring(z.read(part))
        heading = ''
        pi = 0
        body = root.find(W+'body') if part.endswith('document.xml') else root
        # walk all paragraphs in order, noting table context
        for p in body.iter(W+'p'):
            # skip paragraphs nested in another paragraph's textbox duplicates: lxml iter visits nested p too; fine
            pi += 1
            fin, dele = ptext(p)
            if p.find('.//'+W+'p') is not None:  # container paragraph (textbox) -> its children are visited separately
                fin2 = ''.join((e.text or '') for e in p.iterchildren(W+'r') for e in e.iter(W+'t'))
                fin = fin2
            sty = p.find(W+'pPr/'+W+'pStyle')
            sty = sty.get(W+'val') if sty is not None else ''
            if sty.lower().startswith('heading') and fin.strip():
                heading = fin.strip()[:80]
            tc = next((a for a in p.iterancestors(W+'tc')), None)
            loc = {'part': part.split('/')[-1], 'para': pi, 'style': sty, 'heading': heading}
            rowtext = ''
            if tc is not None:
                tr = tc.getparent(); tbl = tr.getparent()
                alltbl = list(root.iter(W+'tbl'))
                loc['table'] = alltbl.index(tbl) + 1
                loc['row'] = list(tbl.iterchildren(W+'tr')).index(tr) + 1
                loc['cell'] = list(tr.iterchildren(W+'tc')).index(tc) + 1
                rowtext = ' | '.join(ptext(c)[0] for c in tr.iterchildren(W+'tc'))
                # header row text
                first = next(tbl.iterchildren(W+'tr'))
                loc['table_header'] = ' | '.join(ptext(c)[0] for c in first.iterchildren(W+'tc'))[:200]
            if fin.strip():
                yield dict(loc=loc, state='final', text=fin, row=rowtext)
            if dele.strip():
                yield dict(loc=loc, state='deleted', text=dele, row=rowtext)

def xlsx_segments(path):
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=False)
    for ws in wb.worksheets:
        rows = list(ws.iter_rows(values_only=True))
        # header row = first row with >=3 non-empty string cells
        hdr_i = next((i for i, r in enumerate(rows[:15]) if sum(isinstance(v, str) and v.strip() != '' for v in r) >= 3), 0)
        hdr = [str(v) if v is not None else '' for v in rows[hdr_i]] if rows else []
        for ri, r in enumerate(rows):
            rt = ' | '.join('' if v is None else str(v) for v in r)
            bins = BIN.findall(rt)
            for ci, v in enumerate(r):
                if v is None: continue
                s = str(v)
                if not s.strip(): continue
                col = openpyxl.utils.get_column_letter(ci + 1)
                yield dict(loc={'sheet': ws.title, 'cell': f'{col}{ri+1}', 'header': hdr[ci] if ci < len(hdr) and ri > hdr_i else ''},
                           state='value', text=s, row=rt[:600], row_bins=bins)

def pptx_segments(path):
    from pptx import Presentation
    prs = Presentation(path)
    for si, sl in enumerate(prs.slides, 1):
        def walk(shapes):
            for sh in shapes:
                if sh.shape_type == 6:  # group
                    yield from walk(sh.shapes)
                if getattr(sh, 'has_text_frame', False) and sh.has_text_frame:
                    t = sh.text_frame.text
                    if t.strip(): yield dict(loc={'slide': si, 'shape': sh.name}, state='text', text=t, row='')
                if getattr(sh, 'has_table', False) and sh.has_table:
                    for ri, row in enumerate(sh.table.rows, 1):
                        rt = ' | '.join(c.text for c in row.cells)
                        for ci, c in enumerate(row.cells, 1):
                            if c.text.strip():
                                yield dict(loc={'slide': si, 'shape': sh.name, 'row': ri, 'cell': ci}, state='table', text=c.text, row=rt)
                if sh.shape_type == 13:
                    yield dict(loc={'slide': si, 'shape': sh.name, 'picture': getattr(sh, 'image', None) and sh.image.filename}, state='picture', text='', row='')
        yield from walk(sl.shapes)
        if sl.has_notes_slide:
            t = sl.notes_slide.notes_text_frame.text
            if t.strip(): yield dict(loc={'slide': si, 'shape': 'notes'}, state='notes', text=t, row='')

def text_segments(path):
    for i, line in enumerate(path.read_text(errors='replace').splitlines(), 1):
        if line.strip():
            yield dict(loc={'line': i}, state='line', text=line, row='')

def pdf_segments(path):
    out = subprocess.run(['pdftotext', '-layout', str(path), '-'], capture_output=True, text=True).stdout
    for pg, page in enumerate(out.split('\f'), 1):
        for i, line in enumerate(page.splitlines(), 1):
            if line.strip():
                yield dict(loc={'page': pg, 'line': i}, state='line', text=line, row='')

def main(files):
    n = 0
    with open(OUT, 'w') as fh:
        for f in files:
            p = EK/f
            ext = p.suffix.lower()
            fn = {'.docx': docx_segments, '.xlsx': xlsx_segments, '.pptx': pptx_segments, '.md': text_segments,
                  '.html': text_segments, '.tsv': text_segments, '.csv': text_segments, '.pdf': pdf_segments}[ext]
            try:
                for s in fn(p):
                    s['file'] = f
                    fh.write(json.dumps(s) + '\n'); n += 1
            except Exception as e:
                print('ERR', f, e, file=sys.stderr)
    print('segments', n)

if __name__ == '__main__':
    main([l.rstrip('\n') for l in open(sys.argv[1]) if l.strip()])
