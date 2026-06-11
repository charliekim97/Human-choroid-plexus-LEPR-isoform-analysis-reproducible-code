#!/usr/bin/env python3
"""
plot_control_main_stacked.py — control-only main figure, Panel c as a 100%-stacked
composition bar (short-form bottom/grey, long-form top/dark) with individual donor
short-fraction points overlaid. Panels a/b unchanged (paired short vs long).
NEW output; does NOT overwrite Fig_LEPR_control_main.*.

In : results/salmon_lepr_tpm.csv, results/LEPR_expression_GSE137619.csv
Out: figures/Fig_LEPR_control_main_stacked.{pdf,png,svg} + _caption.txt
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
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures", "Fig_LEPR_control_main_stacked")
COH = {"disc": "#4C72B0", "repl": "#DD8452"}
SHORT_C, LONG_C = "#cfcfcf", "#3a3a3a"

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

fig, ax = plt.subplots(1, 3, figsize=(10.2, 3.3), gridspec_kw=dict(wspace=0.5))

def paired(a, s, l, color, title, fold, nd, n):
    for i in range(len(s)): a.plot([0, 1], [s[i], l[i]], color="#cccccc", lw=0.8, zorder=1)
    a.scatter(np.zeros(len(s)), s, s=30, c=color, edgecolor="white", lw=0.4, zorder=3)
    a.scatter(np.ones(len(l)), l, s=30, c=color, edgecolor="white", lw=0.4, zorder=3)
    a.set_yscale("log"); a.set_xlim(-0.45, 1.45); a.set_ylim(0.1, 200)
    a.set_xticks([0, 1]); a.set_xticklabels(["Short LEPR\nisoforms", "Long LEPR\nisoform"])
    a.set_ylabel("LEPR isoform abundance (TPM)")
    a.text(0.5, 0.95, f"{nd}/{n} donors · ≈{fold:.0f}-fold", transform=a.get_xaxis_transform(),
           ha="center", va="top", fontsize=7.8, color="#333")
    a.set_title(title, fontsize=8.5)

paired(ax[0], dcs, dcl, COH["disc"], "GSE228458 control (n=5)", dfold, dnd, dn)
ax[0].text(-0.32, 1.06, "a", transform=ax[0].transAxes, fontsize=12, fontweight="bold")
paired(ax[1], rcs, rcl, COH["repl"], "GSE137619 control (n=6)", rfold, rnd, rn)
ax[1].text(-0.32, 1.06, "b", transform=ax[1].transAxes, fontsize=12, fontweight="bold")

# ===== Panel c: 100%-stacked composition (median) + individual donor short-fraction points =====
c = ax[2]; rng = np.random.default_rng(3); W = 0.56
for xpos, sf, ccol, n in [(0, dsf, COH["disc"], dn), (1, rsf, COH["repl"], rn)]:
    msf = np.median(sf)
    c.bar(xpos, msf, width=W, color=SHORT_C, edgecolor="#222", lw=0.7, zorder=1)
    c.bar(xpos, 100-msf, width=W, bottom=msf, color=LONG_C, edgecolor="#222", lw=0.7, zorder=1)
    x = rng.normal(xpos, 0.05, len(sf))
    c.scatter(x, sf, s=30, c=ccol, edgecolor="white", lw=0.6, zorder=4)   # each donor's short fraction
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
(a) GSE228458 control donors (n=5) and (b) GSE137619 control donors (n=6): short LEPR isoform
versus long LEPR isoform abundance (Salmon TPM, log scale); thin lines connect paired values
within each donor; colour denotes cohort. Short-form abundance exceeds the long isoform in
every control donor of both cohorts ({dnd}/{dn} and {rnd}/{rn}; median {dfold:.0f}- and {rfold:.0f}-fold).
(c) LEPR isoform composition in control donors of each cohort, shown as a 100%-stacked bar of
the median short-form (grey) and long-form (dark) fractions; overlaid points are the short-form
fraction of each individual donor (n={dn} and n={rn}). Median short-form fraction was
{np.median(dsf):.1f}% (GSE228458) and {np.median(rsf):.1f}% (GSE137619); one GSE228458 control
(87.7%) is lower but quality-normal (see QC) and still short-dominant.
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
print(f"ctrl medians short-form: GSE228458 {np.median(dsf):.1f}% / GSE137619 {np.median(rsf):.1f}%")
print("wrote Fig_LEPR_control_main_stacked (pdf/png/svg) + caption (original files untouched)")
