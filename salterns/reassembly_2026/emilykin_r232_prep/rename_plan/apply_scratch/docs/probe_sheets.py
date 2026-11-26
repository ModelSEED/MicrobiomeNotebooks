import openpyxl
from collections import Counter
R="/home/freiburger/Documents/EmilyKin"
B="/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scratch/data/backups/MAG_target_gene_filled.HEAD.xlsx"
sub=openpyxl.load_workbook(f"{R}/MAG_target_gene_filled_UPDATED (1).xlsx")
bh=openpyxl.load_workbook(B)
for sh in sub.sheetnames:
    ws=sub[sh]
    print("==",sh,ws.max_row,ws.max_column)
    for r in range(1,5):
        vals=[ws.cell(r,c).value for c in range(1,min(ws.max_column,40)+1)]
        if any(v is not None for v in vals): print(" r%d"%r, [str(v)[:25] if v is not None else None for v in vals])
for f,sh in [("cladeexcel.xlsx",None),("abundance and coverage (1).xlsx",None)]:
    wb=openpyxl.load_workbook(f"{R}/{f}")
    for s in wb.sheetnames:
        ws=wb[s]; print("==",f,s,ws.max_row,ws.max_column)
        for r in range(1,5):
            vals=[ws.cell(r,c).value for c in range(1,min(ws.max_column,30)+1)]
            if any(v is not None for v in vals): print(" r%d"%r, [str(v)[:25] if v is not None else None for v in vals])
