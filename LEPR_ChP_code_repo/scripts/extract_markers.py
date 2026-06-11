#!/usr/bin/env python3
"""ChP contamination QC — extract per-sample marker-gene TPM + LEPR isoform TPM
from Salmon quant.sf. Gene-level TPM = sum of transcript TPM by gene_name
(== tximport gene 'abundance'; matches GENCODE v44 / GRCh38.p14 used for quant).

Usage:
  extract_markers.py <quant_dir> <sample_map.csv> <tx2gene.tsv> <out.csv>
    quant_dir : dir containing <donor>/quant.sf
    tx2gene   : transcript_id <tab> gene_id <tab> gene_name  (built by sbatch)
Outputs one CSV per cohort. Cohorts are NEVER merged downstream (batch/platform differ).
"""
import sys, os, glob, csv, json

quant_dir, smap_path, tx2gene_path, out = sys.argv[1:5]

# ---- marker panel (gene symbols). FOXJ1 deliberately EXCLUDED (ChP+ependyma both +). ----
MARKERS = {
 "ChP":    ["TTR","FOLR1","PRLR","OTX2","AQP1","ENPP2","KCNJ13","CLIC6","HTR2C","KL"],
 "Neuron": ["RBFOX3","SNAP25","SYP","MAP2","STMN2","SLC17A7","GAD1","GAD2"],
 "Endo":   ["CLDN5","PECAM1","FLT1","CDH5"],
 "Astro":  ["GFAP","AQP4","ALDH1L1"],
 "Micro":  ["AIF1","P2RY12","PTPRC"],
}
ALLG = [g for v in MARKERS.values() for g in v]
LEPR_SHORT = {"ENST00000371060","ENST00000616738"}   # 2-tx LepRa def (acceptor 65,633,152)
LEPR_LONG  = {"ENST00000349533"}                      # LepRb 1165aa (acceptor 65,636,191)

# transcript(no version) -> gene_name
tx2name = {}
for line in open(tx2gene_path):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 3:
        tx2name[p[0].split(".")[0]] = p[2]

smap = {r["donor"]: r["group"] for r in csv.DictReader(open(smap_path))}

rows = []
for q in sorted(glob.glob(os.path.join(quant_dir, "*", "quant.sf"))):
    donor = os.path.basename(os.path.dirname(q))
    grp = smap.get(donor, "NA")
    gtpm = {g: 0.0 for g in ALLG}
    short = long = 0.0
    for r in csv.DictReader(open(q), delimiter="\t"):
        tid = r["Name"].split(".")[0]
        tpm = float(r["TPM"])
        gn = tx2name.get(tid)
        if gn in gtpm: gtpm[gn] += tpm
        if tid in LEPR_SHORT: short += tpm
        if tid in LEPR_LONG:  long += tpm
    pm = ""
    mi = os.path.join(os.path.dirname(q), "aux_info", "meta_info.json")
    if os.path.exists(mi):
        try: pm = json.load(open(mi)).get("percent_mapped", "")
        except Exception: pm = ""
    rows.append((donor, grp, pm, gtpm, short, long))

with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["donor", "group", "percent_mapped"] + ALLG + ["LEPR_short_TPM_2tx", "LEPR_long_TPM"])
    for donor, grp, pm, gtpm, short, long in rows:
        w.writerow([donor, grp, pm] + [f"{gtpm[g]:.4f}" for g in ALLG] + [f"{short:.4f}", f"{long:.4f}"])
print("wrote", out, "n=", len(rows))
