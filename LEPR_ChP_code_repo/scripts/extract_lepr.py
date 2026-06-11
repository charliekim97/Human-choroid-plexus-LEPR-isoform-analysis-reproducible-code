#!/usr/bin/env python3
"""
Extract TOTAL LEPR expression (gene-level) + transcript-level TPM from Salmon
quant.sf outputs, conservatively annotate LepRa-like (short) vs LepRb-like (long),
and write LEPR_expression_FINAL.csv + a group summary.

Run from $WORK (where quant/<donor>/quant.sf and sample_map.csv live).
Reference: GENCODE v44 (GRCh38.p14). LEPR gene ENSG00000116678.
"""
import os, glob, csv, statistics as st

# ---- LEPR transcript models (GENCODE v44 / Ensembl), conservative annotation ----
# protein length drives the short/long call:
#   1165 aa  = LepRb (LONG, intracellular signaling / STAT3)  -> LepRb-like
#   ~896-906 = short cytoplasmic tail isoforms                -> LepRa-like
#   958 aa / non-coding (NCD) / retained                       -> 'other' (NOT counted in a/b)
LEPR_TX = {
 'ENST00000349533':('LEPR-202',1165,'LepRb-like_long'),
 'ENST00000344610':('LEPR-201',906,'LepRa-like_short'),
 'ENST00000371058':('LEPR-203',906,'LepRa-like_short'),
 'ENST00000371060':('LEPR-205',896,'LepRa-like_short'),
 'ENST00000616738':('LEPR-208',896,'LepRa-like_short'),
 'ENST00000371059':('LEPR-204',958,'other'),
 'ENST00000462765':('LEPR-206',None,'other_noncoding'),
 'ENST00000471762':('LEPR-207',None,'other_noncoding'),
}
SHORT={k for k,v in LEPR_TX.items() if v[2]=='LepRa-like_short'}
LONG ={k for k,v in LEPR_TX.items() if v[2]=='LepRb-like_long'}

smap={}
for r in csv.DictReader(open('sample_map.csv')):
    smap[r['donor']]=r['group']

rows=[]
for q in sorted(glob.glob('quant/*/quant.sf')):
    donor=q.split('/')[-2]; grp=smap.get(donor,'NA')
    gene_tpm=gene_cnt=0.0; tx_tpm={}; tx_cnt={}
    for r in csv.DictReader(open(q),delimiter='\t'):
        tid=r['Name'].split('.')[0]
        if tid in LEPR_TX:
            tpm=float(r['TPM']); cnt=float(r['NumReads'])
            gene_tpm+=tpm; gene_cnt+=cnt; tx_tpm[tid]=tpm; tx_cnt[tid]=cnt
    short_tpm=sum(tx_tpm.get(t,0) for t in SHORT)
    long_tpm =sum(tx_tpm.get(t,0) for t in LONG)
    frac_long = long_tpm/(short_tpm+long_tpm) if (short_tpm+long_tpm)>0 else float('nan')
    rows.append(dict(donor=donor,group=grp,gene_tpm=gene_tpm,gene_cnt=gene_cnt,
                     short_tpm=short_tpm,long_tpm=long_tpm,frac_long=frac_long,tx_tpm=tx_tpm))

# ---- write per-donor CSV ----
txcols=list(LEPR_TX.keys())
with open('LEPR_expression_FINAL.csv','w',newline='') as f:
    w=csv.writer(f)
    head=['donor','group','LEPR_gene_TPM','LEPR_gene_NumReads',
          'LepRa_like_short_TPM','LepRb_like_long_TPM','frac_long']
    head+=[f"{LEPR_TX[t][0]}_{t}_TPM" for t in txcols]
    w.writerow(head)
    for r in rows:
        row=[r['donor'],r['group'],f"{r['gene_tpm']:.4f}",f"{r['gene_cnt']:.2f}",
             f"{r['short_tpm']:.4f}",f"{r['long_tpm']:.4f}",
             f"{r['frac_long']:.4f}" if r['frac_long']==r['frac_long'] else 'NA']
        row+=[f"{r['tx_tpm'].get(t,0):.4f}" for t in txcols]
        w.writerow(row)

# ---- group summary + control-vs-DM1 test ----
def grp(col):
    out={}
    for g in ['control','AD','DM1']:
        out[g]=[r[col] for r in rows if r['group']==g]
    return out
print("=== TOTAL LEPR gene TPM by group ===")
gt=grp('gene_tpm')
for g in ['control','AD','DM1']:
    v=gt[g]; print(f"  {g}: n={len(v)} mean={st.mean(v):.3f} median={st.median(v):.3f}")
try:
    from scipy.stats import mannwhitneyu
    u,p=mannwhitneyu(gt['control'],gt['DM1'],alternative='two-sided')
    print(f"  control vs DM1 (gene TPM): Mann-Whitney p={p:.3f}")
    sh=grp('short_tpm')
    u2,p2=mannwhitneyu(sh['control'],sh['DM1'],alternative='two-sided')
    print(f"  control vs DM1 (LepRa-like short TPM): p={p2:.3f}")
except ImportError:
    print("  (pip install scipy for p-values)")
print("\nwrote LEPR_expression_FINAL.csv")
print("KEY QUESTION: compare LEPR_gene_TPM and LepRa_like_short_TPM, control vs DM1.")
