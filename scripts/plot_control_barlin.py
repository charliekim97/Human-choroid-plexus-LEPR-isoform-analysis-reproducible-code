#!/usr/bin/env python3
"""
plot_control_barlin.py — control-only, 2-panel (no panel c), LRP2 style:
normalised BAR graph, LINEAR y-axis ending at 0 and 1 (short-form = 1).
Each donor normalised to its OWN short-form (short=1), so nothing exceeds 1 and the
y-axis runs 0 -> 1. Long isoform appears as the small fraction it is (~3-6% of short).
Bars = mean, error = SEM, individual donor points overlaid.
NEW output; does NOT overwrite Fig_LEPR_control_main*, *_stacked, *_barnorm.

In : results/salmon_lepr_tpm.csv, results/LEPR_expression_GSE137619.csv
Out: figures/Fig_LEPR_control_barlin.{pdf,png,svg} + _caption.txt
"""
import csv, os, numpy as np
from scipy.stats import wilcoxon
import matplotlib as mpl; mpl.use("Agg"); import matplotlib.pyplot as plt
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8, "axes.edgecolor": "#222",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures", "Fig_LEPR_control_barlin")
COH = {"disc": "#4C72B0", "repl": "#DD8452"}

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

fig, ax = plt.subplots(1, 2, figsize=(6.6, 3.3), gridspec_kw=dict(wspace=0.45))

def barlin(a, s, l, color, title, nd, n):
    sn, ln = s/s, l/s                       # per-donor: short=1, long=long/short
    means = [np.mean(sn), np.mean(ln)]      # short=1.0 ; long=mean fraction of short
    sems  = [0.0, ln.std(ddof=1)/np.sqrt(len(ln))]
    a.bar([0, 1], means, width=0.6, color=color, alpha=0.30, edgecolor=color, lw=1.1, zorder=1)
    a.errorbar([0, 1], means, yerr=sems, fmt="none", ecolor="#222", elinewidth=0.9, capsize=3, zorder=4)
    rng = np.random.default_rng(1)
    a.scatter(rng.normal(0, 0.07, len(sn)), sn, s=22, c=color, edgecolor="white", lw=0.4, zorder=5)
    a.scatter(rng.normal(1, 0.07, len(ln)), ln, s=22, c=color, edgecolor="white", lw=0.4, zorder=5)
    a.set_ylim(0, 1.05); a.set_yticks([0, 0.5, 1.0]); a.set_xlim(-0.5, 1.5)
    a.set_xticks([0, 1]); a.set_xticklabels(["Short LEPR\nisoforms", "Long LEPR\nisoform"])
    a.set_ylabel("LEPR abundance\n(relative to short-form)")
    a.text(0.97, 0.93, f"{nd}/{n} donors", transform=a.transAxes,
           ha="right", va="top", fontsize=7.8, color="#333")
    a.set_title(title, fontsize=8.5)
    return np.mean(ln), sems[1]

dlong_m, dlong_se = barlin(ax[0], dcs, dcl, COH["disc"], "GSE228458 control (n=5)", dnd, dn)
ax[0].text(-0.30, 1.04, "a", transform=ax[0].transAxes, fontsize=12, fontweight="bold")
rlong_m, rlong_se = barlin(ax[1], rcs, rcl, COH["repl"], "GSE137619 control (n=6)", rnd, rn)
ax[1].text(-0.30, 1.04, "b", transform=ax[1].transAxes, fontsize=12, fontweight="bold")

fig.suptitle("Short-form LEPR predominates in control human choroid plexus (two independent cohorts)",
             fontsize=10, fontweight="bold", y=1.04)
for ext in ("pdf", "png", "svg"): fig.savefig(f"{FIG}.{ext}", dpi=400, bbox_inches="tight")

cap = f"""Figure. Short-form LEPR predominates in control human choroid plexus, replicated across two independent cohorts.
(a) GSE228458 control donors (n=5) and (b) GSE137619 control donors (n=6): LEPR short- and
long-isoform abundance, each donor normalised to its own short-form level (short-form = 1),
linear scale; bars, group mean; error bars, SEM (the short-form bar is the per-donor reference
and has no spread); points, individual donors. The long isoform is a mean of {dlong_m*100:.1f}%
(GSE228458) and {rlong_m*100:.1f}% (GSE137619) of the short-form level; short-form exceeds the
long isoform in every control donor ({dnd}/{dn} and {rnd}/{rn}; median {dfold:.0f}- and {rfold:.0f}-fold).
Statistics: paired Wilcoxon signed-rank (short vs long isoform abundance, donor as unit)
P = {dp:.3f} (GSE228458 control; for n=5 the two-sided exact minimum is 0.0625, so n=5 cannot
reach P<0.05 even with perfect separation) and P = {rp:.3f} (GSE137619 control); pooling both
control sets (n={cn}) gives {cnd}/{cn} short>long, P = {cp:.1g}. Cohorts are independent and not
pooled. "Control" denotes the unaffected comparison group of each source study (GSE228458 =
controls of an Alzheimer disease / myotonic dystrophy type 1 study; GSE137619 = controls of a
progressive multiple sclerosis study; definitions in Methods). Disease-group donors (AD, DM1,
MS) were also analysed and showed concordant short-form predominance (Supplementary).
Short-form isoforms correspond to LepRa-type transcripts; the long isoform to canonical LepRb.
Salmon v1.10.0; GENCODE v44 (GRCh38.p14) decoy-aware; full libraries; paired-end ~150 bp;
donor as the unit; exploratory."""
open(f"{FIG}_caption.txt", "w").write(cap)
print(f"long as % of short: GSE228458 mean {dlong_m*100:.1f}% ; GSE137619 mean {rlong_m*100:.1f}%")
print("wrote Fig_LEPR_control_barlin (pdf/png/svg) + caption (originals untouched)")
