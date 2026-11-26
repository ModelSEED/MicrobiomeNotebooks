"""Bin-keyed cell plan for tabular documents (xlsx sheets, Table S6 tsv, Table S7 csv).
For each row keyed by a MAG bin: MAG-ID cells (era detection + R232 ID), GTDB rank / lineage cells (R220 -> R232),
and GTDB-version / skani-reference cells (flag for recompute). Output: table_cells_plan.tsv + table_sheet_summary.tsv"""
import csv, json, re, sys
from collections import Counter, defaultdict
import openpyxl
EK='/home/freiburger/Documents/EmilyKin/'
RES='/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/results/'
H=json.load(open('id_history.json'))
ERAS={'HEAD':H['head']}
for c in ['ffe668c','6a04588','924fa18']:
    ERAS[c]={m:h.get(c,'') for m,h in H['hist'].items()}
NEW=H['new']
R220={r['user_genome']:r['classification'] for r in csv.DictReader(open(EK+'gtdbtk.bac120.summary.tsv'),delimiter='\t')}
R232={r['user_genome']:r['classification'] for r in csv.DictReader(open(RES+'gtdbtk.r232.summary.tsv'),delimiter='\t')}
RK=['domain','phylum','class','order','family','genus','species']
def ranks(lin): return {RK[i]:x[3:] for i,x in enumerate(lin.split(';'))}
BIN=re.compile(r'^(?:CAN_\d_bin\.\d+|coasm_bin\.\d+)$')
def norm(h): return (h or '').strip().lower()
IDH=lambda h: ('iteration id' in norm(h) or norm(h) in ('iterative id','iterative_id','iterativeid','mag_iterativeid'))
LINH=lambda h: norm(h) in ('taxonomy','gtdb taxonomy','mag_gtdb_taxonomy','taxonomy_full')
RANKH=lambda h: norm(h) in RK
REFH=lambda h: norm(h) in ('reference database','skani_best_taxonomy','closest reference genome name','skani_best_ref','closest reference accession',
                         'ani to closest reference (%)','alignment fraction / af','species-level call','novelty interpretation','novelty_interpretation','ani_class','skani_ani','skani_af_query','skani_af_ref','notes')
cells=[]; summ=[]
def handle_table(fname, sheet, rows, cellname):
    # find header row: first row containing a bin in some col below it -> header = previous non-empty row with strings
    hdr_i=None; bincol=None
    for i,r in enumerate(rows):
        for j,v in enumerate(r):
            if isinstance(v,str) and BIN.match(v.strip()):
                bincol=j; hdr_i=i-1; break
        if bincol is not None: break
    if bincol is None: return
    k=hdr_i
    while k>=0 and not (bincol<len(rows[k]) and isinstance(rows[k][bincol],str) and norm(rows[k][bincol]) in ('mag','mags','mag bin name','bin')): k-=1
    if k>=0: hdr_i=k
    else:
        while hdr_i>=0 and not any(isinstance(v,str) and v.strip() for v in rows[hdr_i]): hdr_i-=1
    hdr=[str(v) if v is not None else '' for v in rows[hdr_i]] if hdr_i>=0 else ['']*len(rows[0])
    era=Counter(); nchg=Counter()
    for i in range(hdr_i+1,len(rows)):
        r=rows[i]
        if bincol>=len(r): continue
        b=r[bincol]
        if not (isinstance(b,str) and BIN.match(b.strip())): continue
        b=b.strip()
        if b not in R232: 
            era['bin_not_in_291']+=1; continue
        o=ranks(R220[b]); n=ranks(R232[b])
        for j,h in enumerate(hdr):
            if j!=bincol and RANKH(h) and (j>=len(r) or r[j] is None or str(r[j]).strip()=='') and n[norm(h)]:
                nchg['taxon_rank_gain']+=1
                cells.append(dict(file=fname,location=f'{sheet}!{cellname(i,j)}',bin=b,column=h,current='(blank)',id_era='',
                    proposed=n[norm(h)],change_type='taxon_rank_gain',context='MAG (bin-keyed row, GTDB rank column; blank in R220)'))
        for j,v in enumerate(r):
            if j==bincol or v is None: continue
            h=hdr[j] if j<len(hdr) else ''
            s=str(v)
            if IDH(h):
                eras=[e for e,mp in ERAS.items() if mp.get(b)==s]
                era['+'.join(eras) or 'no_era_match']+=1
                if s!=NEW[b]:
                    kind='mag_id_change' if H['head'][b]!=NEW[b] else 'mag_id_drift_only'
                    nchg[kind]+=1
                    cells.append(dict(file=fname,location=f'{sheet}!{cellname(i,j)}',bin=b,column=h,current=s,
                        id_era='+'.join(eras) or '?',proposed=NEW[b],change_type=kind,context='MAG (bin-keyed row, MAG ID column)'))
            elif RANKH(h):
                rk=norm(h)
                if s.strip()!=n[rk]:
                    typ='taxon_rank_change' if s.strip()==o[rk] else 'taxon_rank_mismatch_vs_r220'
                    nchg[typ]+=1
                    cells.append(dict(file=fname,location=f'{sheet}!{cellname(i,j)}',bin=b,column=h,current=s,id_era='',
                        proposed=n[rk] or '(blank in R232)',change_type=typ,context='MAG (bin-keyed row, GTDB rank column)'))
                elif o[rk]!=n[rk] and not s.strip():
                    pass
            elif LINH(h) and s.startswith('d__'):
                if s.strip()!=R232[b]:
                    typ='lineage_change' if s.strip()==R220[b] else 'lineage_mismatch_vs_r220'
                    nchg[typ]+=1
                    cells.append(dict(file=fname,location=f'{sheet}!{cellname(i,j)}',bin=b,column=h,current=s,id_era='',
                        proposed=R232[b],change_type=typ,context='MAG (bin-keyed row, GTDB lineage column)'))
            elif REFH(h):
                if re.search(r'R220|r220|2\.7\.2|d__|GC[AF]_|\d',s):
                    nchg['skani_or_version_cell']+=1
                    cells.append(dict(file=fname,location=f'{sheet}!{cellname(i,j)}',bin=b,column=h,current=s[:160],id_era='',
                        proposed='RECOMPUTE (R232 skani / GTDB-Tk 2.8.0 provenance)',change_type='skani_or_version_cell',context='MAG (bin-keyed row, R220-derived reference/ANI column)'))
    # also rows without bin: blank rows / notes with GTDB text
    for i,r in enumerate(rows):
        for j,v in enumerate(r):
            if isinstance(v,str) and re.search(r'R220|r220|2\.7\.2|release ?220|Patescibacteria\b',v) and not (i>hdr_i and bincol<len(r) and isinstance(r[bincol],str) and BIN.match(r[bincol].strip())):
                cells.append(dict(file=fname,location=f'{sheet}!{cellname(i,j)}',bin='',column='(header/note)',current=v[:300],id_era='',
                    proposed='edit text (GTDB-Tk 2.8.0 / R232 / Patescibacteriota)',change_type='header_or_note_text',context='table header/note'))
                nchg['header_or_note_text']+=1
    summ.append(dict(file=fname,sheet=sheet,bin_col=hdr[bincol] if bincol<len(hdr) else '',id_eras=json.dumps(dict(era)),changes=json.dumps(dict(nchg))))

files=sys.argv[1:]
for f in files:
    p=EK+f
    if f.endswith('.xlsx'):
        wb=openpyxl.load_workbook(p,read_only=True,data_only=False)
        for ws in wb.worksheets:
            rows=[list(r) for r in ws.iter_rows(values_only=True)]
            handle_table(f,ws.title,rows,lambda i,j:f'{openpyxl.utils.get_column_letter(j+1)}{i+1}')
    else:
        d='\t' if f.endswith('.tsv') else ','
        rows=list(csv.reader(open(p),delimiter=d))
        handle_table(f,'-',rows,lambda i,j:f'row{i+1}:col{j+1}')
with open('table_cells_plan.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=['file','location','bin','column','current','id_era','proposed','change_type','context'],delimiter='\t'); w.writeheader(); w.writerows(cells)
with open('table_sheet_summary.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(summ[0]),delimiter='\t'); w.writeheader(); w.writerows(summ)
for s in summ: print(s['file'][:40],'|',s['sheet'][:28],'|',s['id_eras'],'|',s['changes'])
