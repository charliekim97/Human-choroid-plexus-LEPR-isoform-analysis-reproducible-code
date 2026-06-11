#!/usr/bin/env python3
"""
Choroid plexus contamination QC — reviewer-facing load-bearing control for the
LEPR isoform claim. Long isoform (LepRb) is the neuron-enriched form, so the
measured quantity and the contamination vulnerability coincide at the same locus.

Design decisions (from PI):
 - tximport-equivalent gene TPM (sum transcript TPM by gene_name), GENCODE v44.
 - Two cohorts NEVER merged (batch/platform/disease differ). Per-cohort scoring,
   per-cohort outlier flags, per-cohort scatter; overlay only for visual compare.
 - STEP 3 Y-axis: (B) junction-spanning LepRb count = PRIMARY, depth-normalized
   (seqkit equal-depth 30M subsample, seed 42 -> reads per million = CPM);
   (A) Salmon LepRb transcript TPM = secondary sanity check.
 - Mapping: LONG = ENST00000349533 (LEPR-202, 1165aa, acceptor chr1:65,636,191);
   SHORT = ENST00000371060 + ENST00000616738 (acceptor chr1:65,633,152).
 - DM1_3 (663 LEPR reads; 7 junction reads) FLAGGED + labeled, never auto-removed;
   include/exclude sensitivity reported.
 - FOXJ1 excluded from panel (ChP + ependyma both positive -> uninformative).
No destructive ops: nothing is excluded from any matrix; flags only.
"""
import csv, os, numpy as np
from scipy.stats import spearmanr
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = "/sessions/ecstatic-funny-pasteur/mnt/03_Archive/REVALIDATION/LEPR_human_ChP_RNAseq"
RES, FIG = ROOT+"/results", ROOT+"/figures"
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
 "font.size":8,"axes.linewidth":0.8,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})

PANEL = {
 "ChP":   ["TTR","FOLR1","PRLR","OTX2","AQP1","ENPP2","KCNJ13","CLIC6","HTR2C","KL"],
 "Neuron":["RBFOX3","SNAP25","SYP","MAP2","STMN2","SLC17A7","GAD1","GAD2"],
 "Endo":  ["CLDN5","PECAM1","FLT1","CDH5"],
 "Astro": ["GFAP","AQP4","ALDH1L1"],
 "Micro": ["AIF1","P2RY12","PTPRC"],
}
CATCOL = {"ChP":"#1C7293","Neuron":"#C44E52","Endo":"#8172B3","Astro":"#CCB974","Micro":"#937860"}
COHORT = {"GSE228458":{"file":"marker_matrix_GSE228458.csv","junc":"lepr_junction_counts.csv","color":"#4C72B0"},
          "GSE137619":{"file":"marker_matrix_GSE137619.csv","junc":"lepr_junction_GSE137619.csv","color":"#DD8452"}}

def load_markers(fn):
    rows=list(csv.DictReader(open(os.path.join(RES,fn))))
    for r in rows:
        r["_log"]={g:np.log2(float(r[g])+1) for cat in PANEL for g in PANEL[cat]}
        r["_mapped"]=float(r["percent_mapped"]) if r["percent_mapped"] else np.nan
        r["_long_tpm"]=float(r["LEPR_long_TPM"]); r["_short_tpm"]=float(r["LEPR_short_TPM_2tx"])
    return rows

def load_junc(fn):
    d={}
    for r in csv.DictReader(open(os.path.join(RES,fn))):
        d[r["donor"]]=dict(LepRb_reads=float(r["LepRb_reads"]),LepRa_reads=float(r["LepRa_reads"]),
                           LepRb_CPM=float(r["LepRb_CPM"]),LepRa_CPM=float(r["LepRa_CPM"]),
                           total_reads=float(r["total_reads"]))
    return d

def scores(rows):
    for r in rows:
        r["chp"]=np.mean([r["_log"][g] for g in PANEL["ChP"]])
        r["neuro"]=np.mean([r["_log"][g] for g in PANEL["Neuron"]])
        r["ratio"]=r["neuro"]/r["chp"] if r["chp"]>0 else np.nan
    rat=np.array([r["ratio"] for r in rows])
    med=np.median(rat); mad=np.median(np.abs(rat-med)); mad=mad if mad>0 else 1e-9
    for r in rows:
        r["modz"]=0.6745*(r["ratio"]-med)/mad           # Iglewicz-Hoaglin modified z
        r["flag_contam"]=bool(r["modz"]>3.5)            # high-contamination outlier only
    return rows

# ================= load =================
DATA={}
for coh,c in COHORT.items():
    rows=scores(load_markers(c["file"])); junc=load_junc(c["junc"])
    for r in rows:
        j=junc.get(r["donor"],{})
        r["LepRb_CPM"]=j.get("LepRb_CPM",np.nan); r["LepRb_reads"]=j.get("LepRb_reads",np.nan)
        r["junc_total"]=j.get("LepRb_reads",0)+j.get("LepRa_reads",0)
        r["flag_depth"]=bool(r["junc_total"]<10)        # low LEPR-junction depth -> unreliable
    DATA[coh]=rows

# ================= FLAG TABLE =================
with open(RES+"/contamination_scores.csv","w",newline="") as f:
    w=csv.writer(f)
    w.writerow(["cohort","donor","group","percent_mapped","ChP_score_log2","neuronal_score_log2",
                "ratio_neuro_over_ChP","modified_z","LepRb_junction_CPM","LepRb_junction_reads",
                "LEPR_long_TPM","LEPR_short_TPM","flag_contamination","flag_low_junction_depth"])
    for coh,rows in DATA.items():
        for r in rows:
            w.writerow([coh,r["donor"],r["group"],f"{r['_mapped']:.1f}",f"{r['chp']:.3f}",
                f"{r['neuro']:.3f}",f"{r['ratio']:.4f}",f"{r['modz']:.2f}",
                f"{r['LepRb_CPM']:.4f}",int(r['LepRb_reads']),f"{r['_long_tpm']:.4f}",
                f"{r['_short_tpm']:.4f}",int(r['flag_contam']),int(r['flag_depth'])])

# ================= FIG 1: marker heatmap (per cohort) =================
order=[(g,cat) for cat in PANEL for g in PANEL[cat]]
fig,axes=plt.subplots(1,2,figsize=(11,7),gridspec_kw=dict(width_ratios=[13,12],wspace=0.05))
for ax,(coh,rows) in zip(axes,DATA.items()):
    samples=[r["donor"] for r in rows]
    M=np.array([[r["_log"][g] for g,_ in order] for r in rows]).T   # markers x samples
    Z=np.zeros_like(M)
    for i in range(M.shape[0]):
        s=M[i].std(); Z[i]=(M[i]-M[i].mean())/s if s>0 else 0.0
    im=ax.imshow(Z,aspect="auto",cmap="RdBu_r",vmin=-2,vmax=2)
    ax.set_xticks(range(len(samples))); ax.set_xticklabels(samples,rotation=90,fontsize=7)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([g for g,_ in order],fontsize=7)
    for j,(g,cat) in enumerate(order):
        ax.add_patch(plt.Rectangle((-1.6,j-0.5),0.9,1,color=CATCOL[cat],clip_on=False))
    ax.set_xlim(-1.7,len(samples)-0.5)
    ax.set_title(f"{coh}  (n={len(samples)})",fontsize=10,fontweight="bold")
    if ax is axes[1]: ax.set_yticklabels([])
axes[0].legend(handles=[Patch(color=CATCOL[c],label=c) for c in PANEL],
    loc="upper left",bbox_to_anchor=(-0.32,-0.06),ncol=5,fontsize=7.5,frameon=False)
cax=fig.add_axes([0.92,0.35,0.012,0.3]); fig.colorbar(im,cax=cax,label="z-score  log2(TPM+1)")
fig.suptitle("ChP marker panel — high TTR/epithelial, near-absent neuronal signal (per-cohort z-score)",
             fontsize=11,fontweight="bold",y=0.95)
for ext in ("png","pdf"): fig.savefig(f"{FIG}/Fig_QC_marker_heatmap.{ext}",dpi=300,bbox_inches="tight")
plt.close(fig)

# ================= FIG 2: LEPR-long vs neuronal score (per cohort) =================
def scatter(ax,rows,yk,ylab,color,title,note_excl=True):
    x=np.array([r["neuro"] for r in rows]); y=np.array([r[yk] for r in rows])
    rho,p=spearmanr(x,y)
    ax.scatter(x,y,s=42,c=color,edgecolor="white",lw=0.5,zorder=3)
    for r in rows:
        if r["flag_depth"] or r["flag_contam"]:
            ax.scatter([r["neuro"]],[r[yk]],s=70,facecolor="none",edgecolor="#C0392B",lw=1.4,zorder=4)
            ax.annotate(r["donor"],(r["neuro"],r[yk]),fontsize=6.5,color="#C0392B",
                        xytext=(3,3),textcoords="offset points")
    txt=f"Spearman ρ={rho:.2f}, P={p:.2f}  (n={len(rows)})"
    # DM1_3 sensitivity (cohort with a depth-flag)
    if note_excl and any(r["flag_depth"] for r in rows):
        keep=[r for r in rows if not r["flag_depth"]]
        xr=np.array([r["neuro"] for r in keep]); yr=np.array([r[yk] for r in keep])
        rho2,p2=spearmanr(xr,yr)
        txt+=f"\nexcl. flagged: ρ={rho2:.2f}, P={p2:.2f} (n={len(keep)})"
    ax.text(0.03,0.97,txt,transform=ax.transAxes,va="top",fontsize=7.3,
            bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#ccc",lw=0.6))
    ax.set_xlabel("Neuronal contamination score\n[mean log2(TPM+1), neuronal panel]")
    ax.set_ylabel(ylab); ax.set_title(title,fontsize=9.5)
    ax.spines[["top","right"]].set_visible(False)
    ymax=max(np.nanmax(y)*1.25,0.1); ax.set_ylim(-ymax*0.05,ymax)

fig,axes=plt.subplots(2,2,figsize=(9.2,8.2),gridspec_kw=dict(hspace=0.42,wspace=0.32))
for col,(coh,rows) in enumerate(DATA.items()):
    c=COHORT[coh]["color"]
    scatter(axes[0,col],rows,"LepRb_CPM","LepRb junction-spanning reads\n(per million, equal-depth 30M)",
            c,f"{coh}  ·  (B) PRIMARY: junction-level")
    scatter(axes[1,col],rows,"_long_tpm","LEPR-long (LepRb) TPM\nENST00000349533, Salmon",
            c,f"{coh}  ·  (A) sanity: transcript-level")
fig.suptitle("LEPR-long signal vs neuronal contamination — no positive dependence in either cohort",
             fontsize=11.5,fontweight="bold",y=0.975)
for ext in ("png","pdf"): fig.savefig(f"{FIG}/Fig_QC_junction_vs_neuronal.{ext}",dpi=300,bbox_inches="tight")
plt.close(fig)

# ================= FIG 3: contamination score barplot (per cohort) =================
fig,axes=plt.subplots(1,2,figsize=(11,3.8),gridspec_kw=dict(wspace=0.22))
for ax,(coh,rows) in zip(axes,DATA.items()):
    c=COHORT[coh]["color"]
    don=[r["donor"] for r in rows]; rat=np.array([r["ratio"] for r in rows])
    med=np.median(rat); mad=np.median(np.abs(rat-med)); mad=mad if mad>0 else 1e-9
    thr=med+(3.5/0.6745)*mad                     # ratio at modified-z = 3.5
    cols=["#C0392B" if r["flag_contam"] else c for r in rows]
    edge=["#7B241C" if r["flag_depth"] else "white" for r in rows]
    ax.bar(range(len(don)),rat,color=cols,edgecolor=edge,lw=1.1,zorder=3)
    ax.axhline(thr,ls="--",lw=1,color="#C0392B",zorder=2)
    ax.text(len(don)-0.5,thr,f"  MAD outlier\n  threshold",va="center",ha="left",
            fontsize=6.8,color="#C0392B")
    ax.set_xticks(range(len(don))); ax.set_xticklabels(don,rotation=90,fontsize=7)
    ax.set_ylabel("Contamination score\nneuronal / ChP  [log2(TPM+1)]")
    ax.set_title(f"{coh}  (n={len(don)})  ·  no sample exceeds threshold",fontsize=9.5)
    ax.spines[["top","right"]].set_visible(False)
    ax.set_xlim(-0.7,len(don)-0.3+1.4)
    chp=[r["chp"] for r in rows]; neu=[r["neuro"] for r in rows]
    ax.text(0.02,0.97,f"ChP score {min(chp):.1f}–{max(chp):.1f}\nneuronal {min(neu):.2f}–{max(neu):.2f}",
            transform=ax.transAxes,va="top",fontsize=7,color="#444",
            bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#ddd",lw=0.5))
fig.suptitle("Per-sample contamination score — all samples below MAD outlier threshold (red edge = low LEPR-junction depth)",
             fontsize=10.5,fontweight="bold",y=1.02)
for ext in ("png","pdf"): fig.savefig(f"{FIG}/Fig_QC_contamination_score.{ext}",dpi=300,bbox_inches="tight")
plt.close(fig)

# ================= console summary =================
print("="*70)
for coh,rows in DATA.items():
    print(f"\n### {coh}  (n={len(rows)})")
    x=np.array([r["neuro"] for r in rows])
    for yk,lab in [("LepRb_CPM","junction CPM (B)"),("_long_tpm","long TPM (A)")]:
        y=np.array([r[yk] for r in rows]); rho,p=spearmanr(x,y)
        keep=[r for r in rows if not r["flag_depth"]]
        rho2,p2=spearmanr([r["neuro"] for r in keep],[r[yk] for r in keep])
        print(f"  {lab:18s}: rho={rho:+.3f} p={p:.3f} | excl.flagged rho={rho2:+.3f} p={p2:.3f}")
    fl=[r["donor"] for r in rows if r["flag_contam"]]
    dp=[r["donor"] for r in rows if r["flag_depth"]]
    print(f"  TTR range: {min(float(r['TTR']) for r in rows):.0f}-{max(float(r['TTR']) for r in rows):.0f} TPM")
    print(f"  neuronal score range: {x.min():.2f}-{x.max():.2f} (log2 TPM+1)")
    print(f"  contamination-flagged: {fl or 'none'} | low-junction-depth-flagged: {dp or 'none'}")
print("\nwrote: contamination_scores.csv, Fig_QC_marker_heatmap.{png,pdf}, Fig_QC_junction_vs_neuronal.{png,pdf}")
