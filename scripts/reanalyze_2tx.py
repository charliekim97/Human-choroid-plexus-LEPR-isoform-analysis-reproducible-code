#!/usr/bin/env python3
"""Re-analysis with narrowed LepRa definition (Dr. Yang's request):
   LepRa(short) = ENST00000371060 + ENST00000616738  (junction-matched LEPR-205/208)
   LepRb(long)  = ENST00000349533
No re-run: recomputed from existing per-transcript Salmon TPM."""
import csv, os, numpy as np
from scipy.stats import wilcoxon
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
ROOT="/sessions/ecstatic-funny-pasteur/mnt/03_Archive/REVALIDATION/LEPR_human_ChP_RNAseq"
RES=ROOT+"/results"; FIGOUT=ROOT+"/figures/Fig_LEPR_control_barlin_2tx"
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
 "font.size":8,"axes.linewidth":0.8,"axes.edgecolor":"#222","xtick.major.width":0.8,"ytick.major.width":0.8,
 "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
COH={"disc":"#4C72B0","repl":"#DD8452"}

def load(fn):
    out=[]
    for r in csv.DictReader(open(os.path.join(RES,fn))):
        s=float(r["LEPR-205_ENST00000371060_TPM"])+float(r["LEPR-208_ENST00000616738_TPM"])
        l=float(r["LEPR-202_ENST00000349533_TPM"])
        out.append(dict(donor=r["donor"],group=r["group"],short=s,long=l))
    return out
d228=load("salmon_lepr_tpm.csv"); d137=load("LEPR_expression_GSE137619.csv")

# ---- write 2-tx per-donor table ----
with open(RES+"/LEPR_2tx_perdonor.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["cohort","donor","group","Short_TPM_2tx","Long_TPM","ShortForm_pct","Long_over_Short"])
    for coh,dd in [("GSE228458",d228),("GSE137619",d137)]:
        for x in dd:
            sf=100*x["short"]/(x["short"]+x["long"]); w.writerow([coh,x["donor"],x["group"],
                f"{x['short']:.3f}",f"{x['long']:.3f}",f"{sf:.1f}",f"{x['long']/x['short']:.4f}"])

def ctrl(dd): return [x for x in dd if x["group"]=="control"]
dc=ctrl(d228); rc=ctrl(d137)
def arr(c,k): return np.array([x[k] for x in c])
dcs,dcl=arr(dc,"short"),arr(dc,"long"); rcs,rcl=arr(rc,"short"),arr(rc,"long")
def stat(s,l):
    _,p=wilcoxon(s,l,alternative="two-sided"); return int(np.sum(s>l)),len(s),p,np.median(s/l)
dnd,dn,dp,dfold=stat(dcs,dcl); rnd,rn,rp,rfold=stat(rcs,rcl)
ccs,ccl=np.concatenate([dcs,rcs]),np.concatenate([dcl,rcl]); _,cp=wilcoxon(ccs,ccl,alternative="two-sided")

# ---- figure (linear normalised bar, 2-panel) ----
fig,ax=plt.subplots(1,2,figsize=(6.6,3.3),gridspec_kw=dict(wspace=0.45))
def barlin(a,s,l,color,title,nd,n):
    sn,ln=s/s,l/s; means=[1.0,np.mean(ln)]; sem=[0.0,ln.std(ddof=1)/np.sqrt(len(ln))]
    a.bar([0,1],means,width=0.6,color=color,alpha=0.30,edgecolor=color,lw=1.1,zorder=1)
    a.errorbar([0,1],means,yerr=sem,fmt="none",ecolor="#222",elinewidth=0.9,capsize=3,zorder=4)
    rng=np.random.default_rng(1)
    a.scatter(rng.normal(0,0.07,len(sn)),sn,s=22,c=color,edgecolor="white",lw=0.4,zorder=5)
    a.scatter(rng.normal(1,0.07,len(ln)),ln,s=22,c=color,edgecolor="white",lw=0.4,zorder=5)
    a.set_ylim(0,1.05); a.set_yticks([0,0.5,1.0]); a.set_xlim(-0.5,1.5)
    a.set_xticks([0,1]); a.set_xticklabels(["Short LEPR\nisoforms","Long LEPR\nisoform"])
    a.set_ylabel("LEPR abundance\n(relative to short-form)")
    a.text(0.97,0.93,f"{nd}/{n} donors",transform=a.transAxes,ha="right",va="top",fontsize=7.8,color="#333")
    a.set_title(title,fontsize=8.5); return np.mean(ln)
dm=barlin(ax[0],dcs,dcl,COH["disc"],"GSE228458 control (n=5)",dnd,dn)
rm=barlin(ax[1],rcs,rcl,COH["repl"],"GSE137619 control (n=6)",rnd,rn)
fig.suptitle("Short-form LEPR predominates in control human choroid plexus (LepRa = ENST00000371060 + ENST00000616738)",
             fontsize=9.5,fontweight="bold",y=1.04)
for ext in ("pdf","png","svg"): fig.savefig(f"{FIGOUT}.{ext}",dpi=400,bbox_inches="tight")
print(f"GSE228458 ctrl: {dnd}/{dn} P={dp:.4f} fold≈{dfold:.0f} long%≈{dm*100:.1f}")
print(f"GSE137619 ctrl: {rnd}/{rn} P={rp:.4f} fold≈{rfold:.0f} long%≈{rm*100:.1f}")
print(f"combined: {int(np.sum(ccs>ccl))}/{len(ccs)} P={cp:.4g}")
print("wrote Fig_LEPR_control_barlin_2tx.{pdf,png,svg} + results/LEPR_2tx_perdonor.csv")
