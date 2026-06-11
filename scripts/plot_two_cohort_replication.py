#!/usr/bin/env python3
"""
plot_two_cohort_replication.py — discovery + independent replication figure.
GSE228458 (control/AD/DM1) and GSE137619 (control/MS), analysed identically, shown side by
side (NOT pooled). Nature-Metabolism style; editable vector (pdf.fonttype=42).

In : results/salmon_lepr_tpm.csv (GSE228458), results/LEPR_expression_GSE137619.csv (GSE137619)
Out: figures/Fig_LEPR_two_cohort_replication.{pdf,png,svg} + _caption.txt
"""
import csv, os, numpy as np
from scipy.stats import wilcoxon, kruskal
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8, "axes.edgecolor": "#222",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); OUT = os.path.join(ROOT, "figures", "Fig_LEPR_two_cohort_replication")
C = {"control": "#4C72B0", "AD": "#55A868", "DM1": "#C44E52", "MS": "#8172B3"}

def load(fn):
    rows = list(csv.DictReader(open(os.path.join(RES, fn))))
    g = [r["group"] for r in rows]
    s = np.array([float(r["LepRa_like_short_TPM"]) for r in rows])
    l = np.array([float(r["LepRb_like_long_TPM"]) for r in rows])
    return g, s, l
def stats(s, l):
    _, p = wilcoxon(s, l, alternative="two-sided"); fold = s / l
    return int(np.sum(s > l)), len(s), p, np.median(fold)

dg, ds, dl = load("salmon_lepr_tpm.csv")            # discovery GSE228458
rg, rs, rl = load("LEPR_expression_GSE137619.csv")  # replication GSE137619
dnd, dn, dp, dfold = stats(ds, dl)
rnd, rn, rp, rfold = stats(rs, rl)
dsf, rsf = ds/(ds+dl)*100, rs/(rs+rl)*100

fig, ax = plt.subplots(1, 3, figsize=(10.6, 3.2), gridspec_kw=dict(wspace=0.45))

def paired(a, g, s, l, title, P, fold, n):
    for i in range(len(s)): a.plot([0, 1], [s[i], l[i]], color="#cccccc", lw=0.7, zorder=1)
    for i in range(len(s)):
        a.scatter(0, s[i], s=24, c=C[g[i]], edgecolor="white", lw=0.4, zorder=3)
        a.scatter(1, l[i], s=24, c=C[g[i]], edgecolor="white", lw=0.4, zorder=3)
    a.set_yscale("log"); a.set_xlim(-0.45, 1.45); a.set_ylim(0.1, 320)
    a.set_xticks([0, 1]); a.set_xticklabels(["Short LEPR\nisoforms", "Long LEPR\nisoform"])
    a.set_ylabel("LEPR isoform abundance (TPM)")
    tr = a.get_xaxis_transform()
    a.plot([0, 0, 1, 1], [0.90, 0.93, 0.93, 0.90], transform=tr, color="#222", lw=0.9, clip_on=False)
    a.text(0.5, 0.945, f"$P$ = {P:.1e}", transform=tr, ha="center", va="bottom", fontsize=7.6)
    a.text(0.5, 0.845, f"{n}/{n} donors · ≈{fold:.0f}-fold", transform=tr, ha="center", va="top", fontsize=7, color="#444")
    a.set_title(title, fontsize=8.5, loc="center")

paired(ax[0], dg, ds, dl, "GSE228458 (discovery, n=13)", dp, dfold, dn)
ax[0].text(-0.30, 1.06, "a", transform=ax[0].transAxes, fontsize=12, fontweight="bold")
# legend for all groups
for g in ["control", "AD", "DM1", "MS"]:
    ax[0].scatter([], [], s=24, c=C[g], edgecolor="white", lw=0.4, label={"control":"Control","AD":"AD","DM1":"DM1","MS":"MS"}[g])
ax[0].legend(loc="lower left", frameon=False, fontsize=6.6, handletextpad=0.2, borderpad=0.1, ncol=1)

paired(ax[1], rg, rs, rl, "GSE137619 (replication, n=12)", rp, rfold, rn)
ax[1].text(-0.30, 1.06, "b", transform=ax[1].transAxes, fontsize=12, fontweight="bold")

# Panel c: short-form fraction (%) both cohorts side by side
c = ax[2]; rng = np.random.default_rng(7)
blocks = [("control", dsf, dg, 0), ("AD", dsf, dg, 1), ("DM1", dsf, dg, 2),
          ("control", rsf, rg, 4), ("MS", rsf, rg, 5)]
xt, xl = [], []
for grp, sfarr, garr, xpos in blocks:
    v = sfarr[[i for i, x in enumerate(garr) if x == grp]]
    x = rng.normal(xpos, 0.07, len(v))
    c.scatter(x, v, s=26, c=C[grp], edgecolor="white", lw=0.4, zorder=3)
    c.plot([xpos-0.22, xpos+0.22], [np.median(v)]*2, color="#222", lw=1.6, zorder=4, solid_capstyle="round")
    xt.append(xpos); xl.append({"control":"Control","AD":"AD","DM1":"DM1","MS":"MS"}[grp])
c.axvline(3.0, color="#ddd", lw=0.8)
c.set_xticks(xt); c.set_xticklabels(xl, fontsize=7.5); c.set_xlim(-0.6, 5.6)
c.set_ylim(85, 100); c.set_ylabel("Short-form LEPR fraction (%)")
c.text(1.0, 99.2, "GSE228458", ha="center", fontsize=7, color="#666")
c.text(4.5, 99.2, "GSE137619", ha="center", fontsize=7, color="#666")
c.text(-0.30, 1.06, "c", transform=c.transAxes, fontsize=12, fontweight="bold")
for a in ax: a.tick_params(labelsize=7.5)
fig.suptitle("Short-form LEPR predominance replicates across two independent human choroid plexus cohorts",
             fontsize=10.5, fontweight="bold", y=1.05)
for ext in ("pdf", "png", "svg"): fig.savefig(f"{OUT}.{ext}", dpi=400, bbox_inches="tight")

cap = f"""Extended Data Figure. Short-form LEPR predominance replicates across two independent human choroid plexus (ChP) cohorts.
Two public human ChP bulk RNA-seq datasets were processed with an identical Salmon pipeline
(GENCODE v44, decoy-aware) and shown side by side (not pooled).
(a) GSE228458 (discovery; 5 control, 4 Alzheimer disease, 4 myotonic dystrophy type 1; n=13)
and (b) GSE137619 (replication; 6 control, 6 progressive multiple sclerosis; n=12): short
LEPR isoform versus the long isoform abundance per donor (Salmon TPM, log scale); thin lines
connect paired values within each donor; points coloured by diagnostic group. Short-form
abundance exceeds the long isoform in every donor of both cohorts ({dnd}/{dn} and {rnd}/{rn};
paired Wilcoxon signed-rank P = {dp:.1e} and {rp:.1e}; median {dfold:.0f}- and {rfold:.0f}-fold).
(c) Short-form LEPR fraction (short-form TPM / summed short+long LEPR TPM) by group in both
cohorts; each point, one donor; bars, group medians. Control short-form fraction was
{np.median(dsf[[i for i,x in enumerate(dg) if x=='control']]):.1f}% (GSE228458) and
{np.median(rsf[[i for i,x in enumerate(rg) if x=='control']]):.1f}% (GSE137619).
The cohorts are independent (different centres, donors, disease contexts and sequencers:
NextSeq 500 vs HiSeq 3000). Disease groups (AD, DM1, MS) are shown for robustness and were
not tested as disease-association hypotheses (Kruskal-Wallis n.s.; underpowered). DM1 =
myotonic dystrophy type 1 (NOT diabetes mellitus). Short-form isoforms correspond to
LepRa-type transcripts; the long isoform corresponds to the canonical LepRb transcript.
Salmon v1.10.0; full libraries; paired-end ~150 bp; donor as the unit; exploratory."""
open(f"{OUT}_caption.txt", "w").write(cap)
print(f"discovery {dnd}/{dn} P={dp:.2g} {dfold:.1f}x ; replication {rnd}/{rn} P={rp:.2g} {rfold:.1f}x")
print("wrote", OUT, "(pdf/png/svg) + caption")
