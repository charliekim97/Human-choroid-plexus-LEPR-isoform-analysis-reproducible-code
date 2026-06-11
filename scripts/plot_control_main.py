#!/usr/bin/env python3
"""
plot_control_main.py — MAIN figure: short-form LEPR predominance in NON-NEUROLOGICAL
CONTROL human choroid plexus, replicated across two independent cohorts (controls only).
Disease groups (AD/DM1/MS) are excluded from the main claim but preserved in a supplementary
summary (notes/Supplementary_disease_groups.md) for transparency.

In : results/salmon_lepr_tpm.csv (GSE228458), results/LEPR_expression_GSE137619.csv (GSE137619)
Out: figures/Fig_LEPR_control_main.{pdf,png,svg} + _caption.txt ; notes/Supplementary_disease_groups.md
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
RES = os.path.join(ROOT, "results"); FIG = os.path.join(ROOT, "figures", "Fig_LEPR_control_main")
SUP = os.path.join(ROOT, "notes", "Supplementary_disease_groups.md")
COH = {"disc": "#4C72B0", "repl": "#DD8452"}   # cohort colours (NOT disease)

def load(fn):
    rows = list(csv.DictReader(open(os.path.join(RES, fn))))
    g = [r["group"] for r in rows]
    s = np.array([float(r["LepRa_like_short_TPM"]) for r in rows])
    l = np.array([float(r["LepRb_like_long_TPM"]) for r in rows])
    return g, s, l
def sub(g, s, l, grp):
    idx = [i for i, x in enumerate(g) if x == grp]; return s[idx], l[idx]
def stat(s, l):
    nd, n = int(np.sum(s > l)), len(s)
    _, p = wilcoxon(s, l, alternative="two-sided"); fold = np.median(s / l)
    return nd, n, p, fold

dg, dS, dL = load("salmon_lepr_tpm.csv")
rg, rS, rL = load("LEPR_expression_GSE137619.csv")
dcs, dcl = sub(dg, dS, dL, "control")   # GSE228458 control (n=5)
rcs, rcl = sub(rg, rS, rL, "control")   # GSE137619 control (n=6)
dnd, dn, dp, dfold = stat(dcs, dcl); rnd, rn, rp, rfold = stat(rcs, rcl)
# combined controls (n=11) — primary, reported in caption (avoids over-weighting the n=5 floor)
ccs = np.concatenate([dcs, rcs]); ccl = np.concatenate([dcl, rcl]); cn = len(ccs)
cnd = int(np.sum(ccs > ccl)); _, cp = wilcoxon(ccs, ccl, alternative="two-sided")
dsf = dcs/(dcs+dcl)*100; rsf = rcs/(rcs+rcl)*100

fig, ax = plt.subplots(1, 3, figsize=(10.0, 3.2), gridspec_kw=dict(wspace=0.5))

def paired(a, s, l, color, title, fold, nd, n):
    for i in range(len(s)): a.plot([0, 1], [s[i], l[i]], color="#cccccc", lw=0.8, zorder=1)
    a.scatter(np.zeros(len(s)), s, s=30, c=color, edgecolor="white", lw=0.4, zorder=3)
    a.scatter(np.ones(len(l)),  l, s=30, c=color, edgecolor="white", lw=0.4, zorder=3)
    a.set_yscale("log"); a.set_xlim(-0.45, 1.45); a.set_ylim(0.1, 200)
    a.set_xticks([0, 1]); a.set_xticklabels(["Short LEPR\nisoforms", "Long LEPR\nisoform"])
    a.set_ylabel("LEPR isoform abundance (TPM)")
    # core message = consistency + effect size (P values in caption; n=5 cannot reach P<0.05)
    a.text(0.5, 0.95, f"{nd}/{n} donors · ≈{fold:.0f}-fold", transform=a.get_xaxis_transform(),
           ha="center", va="top", fontsize=7.8, color="#333")
    a.set_title(title, fontsize=8.5)

paired(ax[0], dcs, dcl, COH["disc"], "GSE228458 control (n=5)", dfold, dnd, dn)
ax[0].text(-0.32, 1.06, "a", transform=ax[0].transAxes, fontsize=12, fontweight="bold")
paired(ax[1], rcs, rcl, COH["repl"], "GSE137619 control (n=6)", rfold, rnd, rn)
ax[1].text(-0.32, 1.06, "b", transform=ax[1].transAxes, fontsize=12, fontweight="bold")

# Panel c: short-form fraction (%) — two cohorts' controls side by side
c = ax[2]; rng = np.random.default_rng(7)
for xpos, v, col in [(0, dsf, COH["disc"]), (1, rsf, COH["repl"])]:
    x = rng.normal(xpos, 0.06, len(v))
    c.scatter(x, v, s=32, c=col, edgecolor="white", lw=0.4, zorder=3)
    c.plot([xpos-0.2, xpos+0.2], [np.median(v)]*2, color="#222", lw=1.8, zorder=4, solid_capstyle="round")
c.axvline(0.5, color="#ddd", lw=0.8)
c.set_xticks([0, 1]); c.set_xticklabels(["GSE228458\ncontrol", "GSE137619\ncontrol"], fontsize=7.5)
c.set_xlim(-0.5, 1.5); c.set_ylim(85, 100); c.set_ylabel("Short-form LEPR fraction (%)")
c.set_title("Controls, both cohorts", fontsize=8.5)
c.text(-0.32, 1.06, "c", transform=c.transAxes, fontsize=12, fontweight="bold")
for a in ax: a.tick_params(labelsize=7.5)
fig.suptitle("Short-form LEPR predominates in control human choroid plexus (two independent cohorts)",
             fontsize=10.5, fontweight="bold", y=1.05)
for ext in ("pdf", "png", "svg"): fig.savefig(f"{FIG}.{ext}", dpi=400, bbox_inches="tight")

cap = f"""Figure. Short-form LEPR predominates in control human choroid plexus, replicated across two independent cohorts.
(a) GSE228458 control donors (n=5) and (b) GSE137619 control donors (n=6): short LEPR isoform
versus long LEPR isoform abundance (Salmon TPM, log scale); thin lines connect paired values
within each donor; colour denotes cohort. Short-form abundance exceeds the long isoform in
every control donor of both cohorts ({dnd}/{dn} and {rnd}/{rn}; median {dfold:.0f}- and {rfold:.0f}-fold).
(c) Short-form LEPR fraction (short-form TPM / summed short+long LEPR TPM) in control donors of
each cohort; each point, one donor; bars, group medians; medians {np.median(dsf):.1f}% and {np.median(rsf):.1f}%.
Statistics: paired Wilcoxon signed-rank (short vs long isoform abundance, donor as unit)
P = {dp:.3f} (GSE228458 control; for n=5 the two-sided exact minimum is 0.0625, so n=5 cannot
reach P<0.05 even with perfect separation) and P = {rp:.3f} (GSE137619 control); pooling both
control sets (n={cn}) gives {cnd}/{cn} short>long, P = {cp:.1g}. Cohorts are independent
(different centres, donors and sequencers) and are not pooled for the per-cohort panels.
"Control" denotes the unaffected comparison group of each source study (GSE228458 = controls of
an Alzheimer disease / myotonic dystrophy type 1 study; GSE137619 = controls of a progressive
multiple sclerosis study; precise definitions in Methods). Disease-group donors (AD, DM1, MS)
were also analysed and showed concordant short-form predominance (Supplementary). Short-form
isoforms correspond to LepRa-type transcripts; the long isoform to the canonical LepRb
transcript. Salmon v1.10.0; GENCODE v44 (GRCh38.p14) decoy-aware; full libraries; paired-end
~150 bp; donor as the unit; exploratory."""
open(f"{FIG}_caption.txt", "w").write(cap)

# ---- Supplementary: disease groups preserved (analysed, concordant) ----
def grp_line(name, g, s, l, key):
    ss, ll = sub(g, s, l, key)
    if len(ss) == 0: return None
    nd, n, p, fold = stat(ss, ll); sf = ss/(ss+ll)*100
    return f"| {name} | {n} | {nd}/{n} | {p:.3f} | {fold:.0f} | {np.median(sf):.1f}% |"
sup = ["# Supplementary — disease-group donors (analysed; excluded from main claim)\n",
       "All disease-group donors were processed by the identical pipeline and show the same",
       "direction (short-form LEPR > long isoform in every donor). They are excluded from the",
       "MAIN figure because the claim concerns normal (non-neurological control) choroid plexus;",
       "they are retained here for transparency (public data; reviewers may ask).\n",
       "| Group | n | short>long | paired Wilcoxon P | median fold | median short-form % |",
       "|---|---|---|---|---|---|"]
for nm, g, s, l, k in [("GSE228458 control", dg, dS, dL, "control"),
                       ("GSE228458 AD", dg, dS, dL, "AD"),
                       ("GSE228458 DM1 (myotonic dystrophy)", dg, dS, dL, "DM1"),
                       ("GSE137619 control", rg, rS, rL, "control"),
                       ("GSE137619 MS", rg, rS, rL, "MS")]:
    ln = grp_line(nm, g, s, l, k)
    if ln: sup.append(ln)
sup += ["\nNote: small per-group n (4-6) limits power; for n=4 and n=5 the two-sided Wilcoxon",
        "exact floors are P=0.125 and P=0.0625 respectively. Direction is concordant throughout.",
        "Per-donor values: results/salmon_lepr_tpm.csv (GSE228458), results/LEPR_expression_GSE137619.csv (GSE137619)."]
open(SUP, "w").write("\n".join(sup) + "\n")
print(f"discovery ctrl {dnd}/{dn} P={dp:.3f} {dfold:.1f}x ; replication ctrl {rnd}/{rn} P={rp:.3f} {rfold:.1f}x")
print("wrote main control figure + caption + Supplementary_disease_groups.md")
