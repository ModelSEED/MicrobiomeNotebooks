#!/usr/bin/env python
"""Measured methane attributable to reconstructed MAGs versus organisms seen only in reads (2026-11-26).

For the six unrestored cores and the equal-growth attribution of pathway_summary.py (attribution.tsv), each guild's
methane is split by the fraction of its read-based biomass that its MAGs carry: MAG biomass / read-based biomass
(specific_rates.tsv). Methanogens: the single methanogen species reconstructed, Methanolobus.1, against all
methanogen mcrA (ANME-related reads excluded, as in pathway_summary.py); MPn consumers: their MAG species against phnJ.
The rest is "reads only": methanogens and MPn consumers not captured by any final MAG (Methanohalophilus,
Methanogenium, Methanococcoides and others, plus Methanolobus beyond its MAG's measured abundance).
Also: the genus split of the methanogen mcrA reads with their placement (kept-MAG contig, other bin, unbinned,
unassembled). Writes attribution_mag_vs_reads.tsv and mcrA_genus_placement.tsv.
Run from salterns/reassembly_2026: ~/Documents/py_venv/bin/python modelling_2026/synthesis/attribution_mag_vs_reads.py
"""
import os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import pandas as pd
src=open("modelling_2026/synthesis/pathway_summary.py").read()
ns={}; exec(src[src.index("METHYLO ="):src.index("mtax = mtax.assign")], ns)
sr=pd.read_csv("modelling_2026/simulations/specific_rates.tsv",sep="\t")
B=sr.pivot_table(index="core",columns="scenario",values="biomass_g_m2")
ch4=sr.groupby("core").CH4_measured_umol_m2_d.first()
at=pd.read_csv("modelling_2026/synthesis/attribution.tsv",sep="\t")
mt=pd.read_csv("read_markers_v2/results/mcrA_taxonomy_unrestored_pooled.tsv",sep="\t")
mt["group"]=mt.genus.fillna("").map(ns["mgroup"])
g=mt.groupby("group").reads.sum(); n_meth=g.drop(["anme","unassigned"],errors="ignore").sum()
frac=n_meth/g.sum()+g.get("unassigned",0)/g.sum()
print("mcrA reads by group:",g.to_dict(),"| methanogen reads:",n_meth)
out={}
for scen in ["mpn_0.07 x1","mpn_0.33 x1","mpn_onlyP x1"]:
    a=at[at.scenario==scen].set_index("core"); rows=[]
    pc={"mpn_0.07 x1":"MePn, ratio 0.07","mpn_0.33 x1":"MePn, ratio 0.33","mpn_onlyP x1":"MePn, only P source"}[scen]
    for c in a.index:
        fM=min(B.loc[c,"Methanogen, methylamines"]/(B.loc[c,"Methanogen, methylamines (read-based abundance)"]*frac),1)
        fP=B.loc[c,pc]/B.loc[c,pc+" (read-based abundance)"]
        meth=a.loc[c,["methylotrophic","hydrogenotrophic","acetoclastic"]].sum()
        rows.append(dict(core=c,CH4_umol_m2_d=ch4[c],MPn=a.loc[c,"MPn"],methanogens=meth,fMAG_methanogen=fM,fMAG_MPn=fP,
            Methanolobus_MAG=meth*fM, MPn_MAGs=a.loc[c,"MPn"]*fP, methanogens_reads_only=meth*(1-fM), MPn_reads_only=a.loc[c,"MPn"]*(1-fP)))
    R=pd.DataFrame(rows).set_index("core"); R["all_MAGs"]=R.Methanolobus_MAG+R.MPn_MAGs; R["all_reads_only"]=1-R.all_MAGs
    out[scen]=R
    print("\n==",scen); print(R.round(3).to_string())
    w=R.drop(columns=["CH4_umol_m2_d"]).mul(R.CH4_umol_m2_d,axis=0).sum()/R.CH4_umol_m2_d.sum()
    print("flux-weighted:", w[["all_MAGs","Methanolobus_MAG","MPn_MAGs","all_reads_only","methanogens_reads_only","MPn_reads_only"]].round(3).to_dict())
gg=mt[mt.group.isin(["methylotrophic","hydrogenotrophic","acetoclastic"])].groupby(["genus","group"])[["reads","kept_mag","other_bin","unbinned","unassembled"]].sum().sort_values("reads",ascending=False)
gg["share"]=gg.reads/n_meth
print(); print(gg.round(3).to_string())
pd.concat({k: v for k, v in out.items()}, names=["scenario"]).to_csv("modelling_2026/synthesis/attribution_mag_vs_reads.tsv", sep="\t")
gg.to_csv("modelling_2026/synthesis/mcrA_genus_placement.tsv", sep="\t")
