#!/usr/bin/env python3
"""Ground-truth 291 statistics for the manuscript update. Deterministic; prints a labeled report
organized by manuscript section. Reuses the exact §2.9 phenotype rules on the 291 KO matrix."""
import csv, glob, json, os, statistics as st
R="."
def num(x):
    try: return float(x)
    except: return None
def load_csv(p):
    with open(p) as f: return list(csv.DictReader(f))
def load_tsv(p):
    with open(p) as f: return list(csv.DictReader(f,delimiter="\t"))
P=print

# ---------- §2.4/3.2 quality + genome ----------
q=load_tsv(f"{R}/mag_quality_summary.tsv")
comp=[num(r["Completeness_pct"]) for r in q]; cont=[num(r["Contamination_pct"]) for r in q]
gs=[num(r["Genome_Size_bp"]) for r in q]; n50=[num(r["N50_bp"]) for r in q]; nc=[num(r["Num_Contigs"]) for r in q]
hq=sum(1 for c,x in zip(comp,cont) if c>=90 and x<5)
P("== QUALITY (§2.4/3.2) ==")
P(f"n_MAGs={len(q)}")
P(f"completeness mean={sum(comp)/len(comp):.2f} median={st.median(comp):.2f} min={min(comp):.2f} max={max(comp):.2f}")
P(f"contamination mean={sum(cont)/len(cont):.2f} median={st.median(cont):.2f} min={min(cont):.2f} max={max(cont):.2f}")
P(f"HQ(>=90 & <5)={hq}  >=90comp={sum(1 for c in comp if c>=90)}  [70,90)={sum(1 for c in comp if 70<=c<90)}  [50,70)={sum(1 for c in comp if 50<=c<70)}")
P(f"genome_size Mbp mean={sum(gs)/len(gs)/1e6:.2f} median={st.median(gs)/1e6:.2f} min={min(gs)/1e6:.2f} max={max(gs)/1e6:.2f}")
P(f"contig_N50 mean={sum(n50)/len(n50):,.0f} median={st.median(n50):,.0f}  Num_Contigs median={st.median(nc):.0f}")
# diffcov breakdown
man=load_tsv(f"{R}/mag_coverage_reanalysis/dereplicated_genomes_update_manifest.tsv")
from collections import Counter
P("diffcov status:", dict(Counter(r["status"] for r in man)))

# ---------- §2.5/3.3 taxonomy ----------
g=load_tsv(f"{R}/gtdbtk.bac120.summary.tsv")
RANKS=["d__","p__","c__","o__","f__","g__","s__"]
def rank(cl,pfx):
    for fld in cl.split(";"):
        fld=fld.strip()
        if fld.startswith(pfx): return fld[len(pfx):]
    return ""
phyla=Counter(rank(r["classification"],"p__") for r in g)
P("\n== TAXONOMY (§2.5/3.3) ==")
P(f"n_MAGs={len(g)}  n_phyla={len([p for p in phyla if p])}")
P("phyla:", dict(sorted(phyla.items(), key=lambda x:-x[1])))
for pfx,name in [("o__","order"),("f__","family"),("g__","genus"),("s__","species")]:
    resolved=sum(1 for r in g if rank(r["classification"],pfx))
    P(f"resolved at {name}: {resolved}")
P(f"distinct genera: {len(set(rank(r['classification'],'g__') for r in g if rank(r['classification'],'g__')))}")
P(f"unnamed species (s__ empty): {sum(1 for r in g if not rank(r['classification'],'s__'))}")
mcol=[c for c in g[0] if 'method' in c.lower() or 'note' in c.lower()]
if 'classification_method' in g[0]:
    P("classification_method:", dict(Counter(r.get('classification_method','') for r in g)))

# ---------- §2.6/3.4 functional panel ----------
# eggNOG-only presence matrix (for union-rescue comparison)
egg=load_csv(f"{R}/ko_presence_absence_matrix.csv"); kos=[c for c in egg[0] if c.startswith("K")]
# UNION presence from mag_abundance_summary (KOfam ∪ eggNOG) — canonical for prevalence + phenotypes
summ=load_tsv(f"{R}/mag_abundance_summary.tsv")
def isone(v): return str(v) in ("1","1.0","True")
pa=[{**{k:("1" if isone(r.get(k,"0")) else "0") for k in kos}, "MAG":r["MAG"]} for r in summ]
prev={k:sum(1 for r in pa if r[k]=="1") for k in kos}
egg_pairs=sum(1 for r in egg for k in kos if r[k]=="1")
uni_pairs=sum(prev.values())
P("\n== FUNCTIONAL PANEL (§2.6/3.4) ==")
P(f"n_MAGs(union)={len(pa)}  n_KOs={len(kos)}")
P(f"UNION MAGxKO pairs={uni_pairs}  eggNOG pairs={egg_pairs}  KOfam-only rescue={uni_pairs-egg_pairs}")
P("per-KO prevalence (UNION):", {k:prev[k] for k in sorted(kos)})
# union stats from eggnog (present matrix = eggNOG) vs kofam tblout
# copy matrix hits
cm=load_csv(f"{R}/ko_copy_number_matrix.csv")
hits=load_csv(f"{R}/ko_hits_long.csv")
P(f"ko_hits_long rows={len(hits)}  unique proteins={len(set(r['protein_id'] for r in hits))}  MAGs with hit={len(set(r['MAG'] for r in hits))}")
# pathway completeness
pw=load_csv(f"{R}/ko_pathway_summary.csv")
paths=["Denitrification","Nitrogen_other","PolyP","Phosphate","PHA","Acetate","Propionate_PHV","Glycogen"]
P("pathway completeness (complete / >=1 / mean%):")
for p in paths:
    tot=f"{p}_genes_total"; pres=f"{p}_genes_present"
    n_complete=sum(1 for r in pw if int(r[pres])==int(r[tot]) and int(r[tot])>0)
    n_ge1=sum(1 for r in pw if int(r[pres])>=1)
    meanpct=sum(int(r[pres])/int(r[tot]) for r in pw if int(r[tot])>0)/len(pw)*100
    P(f"  {p}: {n_complete} / {n_ge1} / {meanpct:.1f}%")

# ---------- §2.9/3.5 phenotypes (rules from §2.9, on 291 KO matrix) ----------
def has(r,k): return r.get(k,"0")=="1"
DEN=["K00370","K00371","K00374","K02567","K02568","K15864","K00368","K04561","K02305","K00376"]
UP=[k for k in DEN if k!="K00376"]
def pheno(r):
    core=has(r,"K00937") or has(r,"K22468")
    ptrans=(has(r,"K02040")and has(r,"K02037")and has(r,"K02038")and has(r,"K02036")) or has(r,"K02035")
    preg=has(r,"K07657")or has(r,"K07636")or has(r,"K02039")
    pha=any(has(r,k) for k in ["K00626","K00023","K03821","K05973"])
    phasyn=has(r,"K00626")and has(r,"K00023")and has(r,"K03821")
    glysyn=has(r,"K00703")and has(r,"K00700")
    glysup=has(r,"K00702")or has(r,"K00688")or has(r,"K01214")
    carb=any(has(r,k) for k in ["K00634","K00925","K01895","K01659","K01720","K01908","K01847","K01848"])
    weak=not core
    no3=(has(r,"K00370")and has(r,"K00371")and has(r,"K00374"))or(has(r,"K02567")and has(r,"K02568"))
    no2=has(r,"K15864")or has(r,"K00368")
    nosz=has(r,"K00376")
    anyden=any(has(r,k) for k in DEN); anyup=any(has(r,k) for k in UP)
    complete=no3 and no2 and has(r,"K04561") and nosz
    d=dict(core_polyP_gene=core,phosphate_transport=ptrans,phosphate_regulation_support=preg,
        PHA_module=pha,PHA_synthesis_module=phasyn,glycogen_synthesis_module=glysyn,
        supportive_glycogen_genes=glysup,carbon_activation_module=carb,weak_absent_polyP=weak,
        complete_denitrification=complete,incomplete_denitrification=(anyden and not complete),
        N2O_sink_denitrifying=(nosz and anyup),N2O_sink_non_denitrifying=(nosz and not anyup),
        respiratory_DNRA_potential=has(r,"K03385"),possible_nitrite_ammonification=(has(r,"K00362")or has(r,"K00363")),
        polyP_capable_PAO_unresolved=core,stronger_polyP_phosphate_cycling=(core and ptrans),
        putative_PAO_like=(core and ptrans and pha),
        stronger_putative_PAO_like=(core and ptrans and pha and preg and phasyn),
        GAO_like=(weak and glysyn and carb),
        putative_GAO_like_storage=(weak and glysyn and carb and glysup and pha),
        PHA_capable_non_glycogen_storage=(weak and pha and not glysyn),
        nosZ_present=nosz, anyden=anyden)
    return d
ph=[pheno(r) for r in pa]
P("\n== PHENOTYPES (§2.9/3.5) ==")
for k in ["core_polyP_gene","phosphate_transport","phosphate_regulation_support","PHA_module","PHA_synthesis_module",
          "glycogen_synthesis_module","supportive_glycogen_genes","carbon_activation_module","weak_absent_polyP",
          "complete_denitrification","incomplete_denitrification","N2O_sink_denitrifying","N2O_sink_non_denitrifying",
          "respiratory_DNRA_potential","possible_nitrite_ammonification","polyP_capable_PAO_unresolved",
          "stronger_polyP_phosphate_cycling","putative_PAO_like","stronger_putative_PAO_like",
          "GAO_like","putative_GAO_like_storage","PHA_capable_non_glycogen_storage","nosZ_present","anyden"]:
    P(f"  {k}={sum(1 for d in ph if d[k])}")

# ---------- guilds (Table S7) ----------
gu=load_csv(f"{R}/mag_guild_classification.csv")
gcol=[c for c in gu[0] if c.lower()=="guild"] or [list(gu[0])[-1]]
P("\n== GUILDS (Table S7) ==", dict(Counter(r[gcol[0]] for r in gu)))

# ---------- §2.7/3.6 linkage + §S6 novelty ----------
am=json.load(open(f"{R}/asv_to_mag_mapping.json"))["asv_to_mag"]
P("\n== ASV->MAG LINKAGE (§2.7/3.6) ==")
P("tiers:", dict(Counter(v["confidence_tier"] for v in am.values())))
P(f"with best hit={sum(1 for v in am.values() if v['best_mag'])}  MAGs hit>=1={len(set(v['best_mag'] for v in am.values() if v['best_mag']))}")
nm=load_csv(f"{R}/mag_iterativeID_old_to_new.csv")
P("naming case:", dict(Counter(r.get('assignment_case','') for r in nm)) if 'assignment_case' in nm[0] else "(no case col)")
sk=load_tsv(f"{R}/Table_S6_skani_novelty.tsv")
P("skani ANI_class:", dict(Counter(r["ANI_class"] for r in sk)))

# ---------- 16S extract ----------
n16=0; m16=set()
for line in open(f"{R}/mag_16s.fasta"):
    if line.startswith(">"): n16+=1; m16.add(line[1:].split("__")[0])
P(f"\n== 16S == seqs={n16} MAGs_with_16S={len(m16)}")

# ---------- nosZ per-bin + ppk1 ----------
pb=load_tsv(f"{R}/clade_classify/out/allbins_nosz_clades_bb_per_bin.tsv")
P("\n== nosZ per-bin (§3.8) ==", f"n_MAGs={len(pb)}", "clade_call:", dict(Counter(r["clade_call"] for r in pb)))
ppk=load_tsv(f"{R}/ppk1_classify/ppk1_classification.tsv")
P("ppk1 (§3.9):", [(r['MAG'],r['ppk1_type']) for r in ppk])

# ---------- Bakta totals (291) ----------
tot=Counter()
for d in glob.glob(f"/scratch1/afreiburger/emilykin/processed/mag/bakta/*/*.txt"):
    mag=os.path.basename(os.path.dirname(d))
    if not os.path.exists(f"{R}/dereplicated_genomes/{mag}.fa"): continue
    for line in open(d):
        for key in ["CDSs","tRNAs","tmRNAs","rRNAs","ncRNAs","ncRNA regions","pseudogenes","sORFs","CRISPR arrays"]:
            if line.startswith(key+":"):
                try: tot[key]+=int(line.split(":")[1].strip())
                except: pass
P("\n== BAKTA totals (291) ==", dict(tot))
