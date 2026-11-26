import json, csv, re
EK='/home/freiburger/Documents/EmilyKin/'
RES='/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/results/'
H=json.load(open('id_history.json'))
R220={r['user_genome']:r['classification'] for r in csv.DictReader(open(EK+'gtdbtk.bac120.summary.tsv'),delimiter='\t')}
R232={r['user_genome']:r['classification'] for r in csv.DictReader(open(RES+'gtdbtk.r232.summary.tsv'),delimiter='\t')}
def gcol(lin):
    p=lin.split(';'); g=p[5][3:]
    return g if g else 'f__'+p[4][3:]
def low2(lin):
    p=[x for x in lin.split(';') if len(x)>3]
    return ';'.join(p[-2:])
out=[]
# main _correction docx Table 1
rows={}
for l in open('segments.jsonl'):
    s=json.loads(l)
    if s['file']=='papert3_070626 (2)_correction.docx' and s['loc'].get('table')==1 and s['state']=='final':
        rows.setdefault(s['loc']['row'],{})[s['loc']['cell']]=s['text']
for rn in sorted(rows):
    c=rows[rn]
    b=c.get(1,'').strip()
    if b not in R232: continue
    idc=c.get(2,'').strip(); gc=c.get(3,'').strip()
    newid=H['new'][b]; g220=gcol(R220[b]); g232=gcol(R232[b])
    if idc!=newid:
        out.append(dict(file='papert3_070626 (2)_correction.docx',location=f'Table 1 row {rn} col 2 (iterativeID)',bin=b,current=idc,proposed=newid,
            note='era=HEAD' if idc==H['head'][b] else f'era mismatch (HEAD {H["head"][b]})'))
    if gc!=g232:
        out.append(dict(file='papert3_070626 (2)_correction.docx',location=f'Table 1 row {rn} col 3 (GTDB genus)',bin=b,current=gc,proposed=g232,
            note='matches R220' if gc==g220 else f'differs from R220 rule ({g220})'))
# Emily Oct-21 PDF Table 1 (bins reconstructed from wrapped layout text)
txt=open('/tmp/claude-1000/-home-freiburger-Documents-MicrobiomeNotebooks-salterns/0f6dc778-e8f5-5027-b895-50697230fd54/scratchpad/p3.txt').read().splitlines()
pdfbins=[]
for i,l in enumerate(txt):
    m=re.match(r'\s*(CAN_\d_bi|coasm_bi)\s',l)
    if m and i+1<len(txt):
        n=re.match(r'\s*n\.(\d+)',txt[i+1])
        if n: pdfbins.append((m.group(1)+'n.'+n.group(1), i+1))
for b,ln in pdfbins:
    if b not in R232: print('??',b); continue
    if R220[b]!=R232[b]:
        out.append(dict(file='paper3/papert3_071626.pdf',location=f'Table 1 (Taxonomy_full), text line ~{ln}',bin=b,current=R220[b],proposed=R232[b],
            note='full lineage string; Emily master (no DOCX on disk)'))
print('main docx Table1 rows',len(rows)-1,'; pdf Table1 bins',len(pdfbins))
with open('table1_plan.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0]),delimiter='\t'); w.writeheader(); w.writerows(out)
for o in out: print(o['file'][:12],'|',o['location'],'|',o['bin'],'|',o['current'][-70:],'->',o['proposed'][-70:],'|',o['note'])
