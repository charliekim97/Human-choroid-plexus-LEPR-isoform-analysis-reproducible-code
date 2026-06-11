#!/usr/bin/env python3
"""
plot_control_main_barnorm.py — control-only main figure, supervisor (LRP2) style:
Panels a/b = NORMALISED BAR graphs (abundance relative to mean short-form = 1.0), LOG y so
the long isoform (~0.04) stays visible; bar = mean, error = SEM, individual donor points
overlaid. Panel c = 100%-stacked composition (short grey / long dark) + donor points.
NEW output; originals (Fig_LEPR_control_main*, _stacked*) are NOT overwritten.

In : results/salmon_lepr_tpm.csv, results/LEPR_expression_GSE137619.csv
Out: figures/Fig_LEPR_control_main_barnorm.{pdf,png,svg} + _caption.txt
"""
import csv, os, numpy as np
from scipy.stats import wilcoxon
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8, "axes.edgecolor": "#222",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures", "Fig_LEPR_control_main_barnorm")
COH = {"disc": "#4C72B0", "repl": "#DD8452"}; SHORT_C, LONG_C = "#cfcfcf", "#3a3a3a"

def load(fn):
    rows = list(csv.DictReader(open(os.path.join(RES, fn))))
    g = [r["group"] for r in rows]
    s = np.array([float(r["LepRa_like_short_TPM"]) for r in rows])
    l = np.array([float(r["LepRb_like_long_TPM"]) for r in rows])
    return g, s, l
def sub(g, s, l, grp):
    idx = [i for i, x in enumerate(g) if x == grp]; return s[idx], l[idx]
def stat(s, l):
    _, p = wilcoxon(s, l, alternative="two-sided"); return int(np.sum(s > l)), len(s), p, np.median(s/l)

dg, dS, dL = load("salmon_lepr_tpm.csv"); rg, rS, rL = load("LEPR_expression_GSE137619.csv")
dcs, dcl = sub(dg, dS, dL, "control"); rcs, rcl = sub(rg, rS, rL, "control")
dnd, dn, dp, dfold = stat(dcs, dcl); rnd, rn, rp, rfold = stat(rcs, rcl)
ccs, ccl = np.concatenate([dcs, rcs]), np.concatenate([dcl, rcl]); cn = len(ccs)
cnd = int(np.sum(ccs > ccl)); _, cp = wilcoxon(ccs, ccl, alternative="two-sided")
dsf, rsf = dcs/(dcs+dcl)*100, rcs/(rcs+rcl)*100

fig, ax = plt.subplots(1, 3, figsize=(10.4, 3.3), gridspec_kw=dict(wspace=0.5))

def normbar(a, s, l, color, title, fold, nd, n):
    m = np.mean(s)                       # normalise to mean short-form of this cohort
    sn, ln = s/m, l/m
    means = [np.mean(sn), np.mean(ln)]
    sems = [sn.std(ddof=1)/np.sqrt(len(sn)), ln.std(ddof=1)/np.sqrt(len(ln))]
    a.bar([0, 1], means, width=0.56, color=color, alpha=0.28, edgecolor=color, lw=1.1, zorder=1)
    a.errorbar([0, 1], means, yerr=sems, fmt="none", ecolor="#222", elinewidth=0.9, capsize=3, zorder=4)
    rng = np.random.default_rng(1)
    a.scatter(rng.normal(0, 0.06, len(sn)), sn, s=22, c=color, edgecolor="white", lw=0.4, zorder=5)
    a.scatter(rng.normal(1, 0.06, len(ln)), ln, s=22, c=color, edgecolor="white", lw=0.4, zorder=5)
    a.set_yscale("log"); a.set_ylim(0.008, 5); a.set_xlim(-0.5, 1.5)
    a.set_xticks([0, 1]); a.set_xticklabels(["Short LEPR\nisoforms", "Long LEPR\nisoform"])
    a.set_ylabel("LEPR abundance\n(relative to mean short-form)")
    a.text(0.5, 0.96, f"{nd}/{n} donors · ≈{fold:.0f}-fold", transform=a.get_xaxis_transform(),
           ha="center", va="top", fontsize=7.6, color="#333")
    a.set_title(title, fontsize=8.5)

normbar(ax[0], dcs, dcl, COH["disc"], "GSE228458 control (n=5)", dfold, dnd, dn)
ax[0].text(-0.34, 1.06, "a", transform=ax[0].transAxes, fontsize=12, fontweight="bold")
normbar(ax[1], rcs, rcl, COH["repl"], "GSE137619 control (n=6)", rfold, rnd, rn)
ax[1].text(-0.34, 1.06, "b", transform=ax[1].transAxes, fontsize=12, fontweight="bold")

# Panel c: 100%-stacked composition + donor short-fraction points
c = ax[2]; rng = np.random.default_rng(3); W = 0.56
for xpos, sf, ccol, n in [(0, dsf, COH["disc"], dn), (1, rsf, COH["repl"], rn)]:
    msf = np.median(sf)
    c.bar(xpos, msf, width=W, color=SHORT_C, edgecolor="#222", lw=0.7, zorder=1)
    c.bar(xpos, 100-msf, width=W, bottom=msf, color=LONG_C, edgecolor="#222", lw=0.7, zorder=1)
    c.scatter(rng.normal(xpos, 0.05, len(sf)), sf, s=28, c=ccol, edgecolor="white", lw=0.6, zorder=4)
    c.text(xpos, 102, f"n={n}", ha="center", va="bottom", fontsize=7.5)
c.set_xticks([0, 1]); c.set_xticklabels(["GSE228458\ncontrol", "GSE137619\ncontrol"], fontsize=7.5)
c.set_xlim(-0.6, 1.6); c.set_ylim(0, 100); c.set_ylabel("LEPR isoform composition (%)")
c.text(-0.32, 1.06, "c", transform=c.transAxes, fontsize=12, fontweight="bold")
c.legend(handles=[mpatches.Patch(facecolor=SHORT_C, edgecolor="#222", lw=0.7, label="Short-form"),
                  mpatches.Patch(facecolor=LONG_C, edgecolor="#222", lw=0.7, label="Long-form")],
         loc="upper center", bbox_to_anchor=(0.5, -0.20), frameon=False, fontsize=6.8, ncol=2, handlelength=1.1)

for a in ax: a.tick_params(labelsize=7.5)
fig.suptitle("Short-form LEPR predominates in control human choroid plexus (two independent cohorts)",
             fontsize=10.5, fontweight="bold", y=1.05)
for ext in ("pdf", "png", "svg"): fig.savefig(f"{FIG}.{ext}", dpi=400, bbox_inches="tight")

cap = f"""Figure. Short-form LEPR predominates in control human choroid plexus, replicated across two independent cohorts.
(a) GSE228458 control donors (n=5) and (b) GSE137619 control donors (n=6): LEPR short- and
long-isoform abundance normalised to the cohort mean short-form level (mean short-form = 1),
on a log scale; bars, group mean; error bars, SEM; points, individual donors. The long isoform
is ~3-4% of the short-form level (short-form exceeds long in {dnd}/{dn} and {rnd}/{rn} donors;
median {dfold:.0f}- and {rfold:.0f}-fold). (c) LEPR isoform composition shown as a 100%-stacked
bar of median short-form (grey) and long-form (dark) fractions; overlaid points, each donor's
short-form fraction (medians {np.median(dsf):.1f}% and {np.median(rsf):.1f}%).
Statistics: paired Wilcoxon signed-rank (short vs long isoform abundance, donor as unit)
P = {dp:.3f} (GSE228458 control; for n=5 the two-sided exact minimum is 0.0625, so n=5 cannot
reach P<0.05 even with perfect separation) and P = {rp:.3f} (GSE137619 control); pooling both
control sets (n={cn}) gives {cnd}/{cn} short>long, P = {cp:.1g}. Cohorts are independent and not
pooled for the per-cohort panels. "Control" denotes the unaffected comparison group of each
source study (GSE228458 = controls of an Alzheimer disease / myotonic dystrophy type 1 study;
GSE137619 = controls of a progressive multiple sclerosis study; definitions in Methods).
Disease-group donors (AD, DM1, MS) were also analysed and showed concordant short-form
predominance (Supplementary). Short-form isoforms correspond to LepRa-type transcripts; the
long isoform to canonical LepRb. Salmon v1.10.0; GENCODE v44 (GRCh38.p14) decoy-aware; full
libraries; paired-end ~150 bp; donor as the unit; exploratory."""
open(f"{FIG}_caption.txt", "w").write(cap)
print("wrote Fig_LEPR_control_main_barnorm (pdf/png/svg) + caption (originals untouched)")
