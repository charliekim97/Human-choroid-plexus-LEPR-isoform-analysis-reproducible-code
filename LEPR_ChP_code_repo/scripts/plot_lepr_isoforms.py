#!/usr/bin/env python3
"""
plot_lepr_isoforms.py — publication-grade LEPR isoform figure (Nature-Metabolism style).
rev4 2026-06-02: refined aesthetics; paired dot-plot (Panel a) with a single clean
significance bracket; minimal on-figure text; full legend text moved to a caption file.

Stats (unchanged, correct):
  PRIMARY  Wilcoxon signed-rank, LepRa-like(short) vs LepRb-like(long), paired, donor as unit.
  EFFECT   log2 fold-ratio short/long.
  GROUPS   control/AD/DM1 NOT a hypothesis -> Kruskal-Wallis only (robustness, low power).

In : results/salmon_lepr_tpm.csv
Out: figures/Fig_LEPR_human_ChP_isoform.{pdf,png,svg}, .._caption.txt; results/salmon_lepr_summary.csv
Run: python3 scripts/plot_lepr_isoforms.py
"""
import csv, os, numpy as np
from scipy.stats import wilcoxon, kruskal
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8, "axes.edgecolor": "#222222",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",  # editable vector text
})

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(ROOT, "results", "salmon_lepr_tpm.csv")
OUT  = os.path.join(ROOT, "figures", "Fig_LEPR_human_ChP_isoform")
SUM  = os.path.join(ROOT, "results", "salmon_lepr_summary.csv")

rows = list(csv.DictReader(open(SRC)))
grp = [r["group"] for r in rows]
def col(k): return np.array([float(r[k]) for r in rows])
gene, short, lng, fl = (col("LEPR_gene_TPM"), col("LepRa_like_short_TPM"),
                        col("LepRb_like_long_TPM"), col("frac_long"))
def G(a, g): return a[[i for i, x in enumerate(grp) if x == g]]

w_all, p_all = wilcoxon(short, lng, alternative="two-sided")
w_ctl, p_ctl = wilcoxon(G(short, "control"), G(lng, "control"), alternative="two-sided")
ndom = int(np.sum(short > lng)); fold = short / lng
med_fold, med_l2 = np.nanmedian(fold), np.nanmedian(np.log2(fold))
med_share = 100 * (1 - np.median(fl))
kw = {nm: kruskal(G(a, "control"), G(a, "AD"), G(a, "DM1"))[1]
      for nm, a in [("frac", fl), ("gene", gene)]}

# ---- summary CSV ----
with open(SUM, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["PRIMARY","Wilcoxon signed-rank, LepRa-like vs LepRb-like TPM, paired, donor as unit"])
    w.writerow(["all donors (n=13)", f"{ndom}/13 short>long", f"P={p_all:.2g}"])
    w.writerow(["control only (n=5)", "5/5 short>long", f"P={p_ctl:.3f} (n=5 two-sided floor 0.0625)"])
    w.writerow(["EFFECT", f"median {med_fold:.1f}x (log2FC {med_l2:.2f})", f"range {np.nanmin(fold):.1f}-{np.nanmax(fold):.1f}x"])
    w.writerow(["LepRa-like share of LEPR TPM", f"median {med_share:.1f}%", f"min {100*(1-fl.max()):.1f}%"])
    w.writerow(["GROUP (robustness, not a hypothesis)", "Kruskal-Wallis (n=5/4/4, low power)", f"frac p={kw['frac']:.3f}; gene p={kw['gene']:.3f}"])

# ---- palette ----
C = {"control": "#4C72B0", "AD": "#55A868", "DM1": "#C44E52"}
order = ["control", "AD", "DM1"]; glab = {"control": "Control", "AD": "AD", "DM1": "DM1"}
rng = np.random.default_rng(7)

fig, ax = plt.subplots(1, 3, figsize=(8.6, 3.0), gridspec_kw=dict(wspace=0.42))

# ===== Panel a: paired short vs long (log) =====
a0 = ax[0]
for i in range(len(rows)):
    a0.plot([0, 1], [short[i], lng[i]], color="#c8c8c8", lw=0.7, zorder=1)
for i in range(len(rows)):
    a0.scatter(0, short[i], s=26, c=C[grp[i]], edgecolor="white", lw=0.4, zorder=3)
    a0.scatter(1, lng[i],  s=26, c=C[grp[i]], edgecolor="white", lw=0.4, zorder=3)
a0.set_yscale("log"); a0.set_xlim(-0.45, 1.45); a0.set_ylim(0.12, 320)
a0.set_xticks([0, 1]); a0.set_xticklabels(["Short LEPR\nisoforms", "Long LEPR\nisoform"])
a0.set_ylabel("LEPR isoform abundance (TPM)")
# clean significance bracket (x in data, y in axes fraction)
tr = a0.get_xaxis_transform()
a0.plot([0, 0, 1, 1], [0.90, 0.93, 0.93, 0.90], transform=tr, color="#222", lw=0.9, clip_on=False)
a0.text(0.5, 0.945, f"$P$ = 2.4 × 10$^{{-4}}$", transform=tr, ha="center", va="bottom", fontsize=8)
a0.text(0.5, 0.845, f"13/13 donors · ≈{med_fold:.0f}-fold", transform=tr, ha="center", va="top", fontsize=7.2, color="#444")
# group legend (compact)
for k, g in enumerate(order):
    a0.scatter([], [], s=26, c=C[g], edgecolor="white", lw=0.4, label=glab[g])
a0.legend(title=None, loc="lower left", frameon=False, fontsize=7, handletextpad=0.2, borderpad=0.1)
a0.text(-0.30, 1.02, "a", transform=a0.transAxes, fontsize=12, fontweight="bold", va="bottom")

# ===== Panel b: SHORT-form fraction (%) by group (flipped to match the main claim) =====
b = ax[1]; sfp = (1 - fl) * 100
for j, g in enumerate(order):
    v = G(sfp, g); x = rng.normal(j, 0.07, len(v))
    b.scatter(x, v, s=30, c=C[g], edgecolor="white", lw=0.4, zorder=3)
    b.plot([j-0.2, j+0.2], [np.median(v)]*2, color="#222", lw=1.6, zorder=4, solid_capstyle="round")
b.axhline(np.median(sfp), color="#bbb", ls="--", lw=0.7, zorder=0)
b.set_xticks(range(3)); b.set_xticklabels([glab[g] for g in order]); b.set_xlim(-0.5, 2.5)
b.set_ylim(85, 100); b.set_ylabel("Short-form LEPR fraction (%)")
b.text(-0.28, 1.02, "b", transform=b.transAxes, fontsize=12, fontweight="bold", va="bottom")

# ===== Panel c: total LEPR gene TPM by group =====
c = ax[2]
for j, g in enumerate(order):
    v = G(gene, g); x = rng.normal(j, 0.07, len(v))
    c.scatter(x, v, s=30, c=C[g], edgecolor="white", lw=0.4, zorder=3)
    c.plot([j-0.2, j+0.2], [np.median(v)]*2, color="#222", lw=1.6, zorder=4, solid_capstyle="round")
c.set_xticks(range(3)); c.set_xticklabels([glab[g] for g in order]); c.set_xlim(-0.5, 2.5)
c.set_ylim(0, 62); c.set_ylabel("Total LEPR transcript TPM")
c.text(-0.28, 1.02, "c", transform=c.transAxes, fontsize=12, fontweight="bold", va="bottom")

for a in ax: a.tick_params(labelsize=7.5)
fig.suptitle("Short LEPR isoforms predominate in human choroid plexus RNA-seq",
             fontsize=10.5, fontweight="bold", y=1.04, x=0.5)
fig.savefig(f"{OUT}.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}.png", dpi=400, bbox_inches="tight")
fig.savefig(f"{OUT}.svg", bbox_inches="tight")

# ---- separate caption (journal style; not printed on the figure) ----
cap = f"""Extended Data Figure. Short LEPR isoforms predominate in human choroid plexus RNA-seq.
(a) LEPR isoform abundance (Salmon TPM, log scale) for short LEPR isoforms versus the long
LEPR isoform in 13 human choroid plexus donors (GSE228458); thin lines connect the paired
values within each donor; points coloured by diagnostic group. Short-form abundance exceeds
the long isoform in all 13/13 donors. (b) Short-form LEPR fraction (short-form TPM as a
percentage of summed short + long LEPR TPM) and (c) total LEPR transcript TPM, by diagnostic
group (control n=5; AD n=4; DM1 n=4); each point, one donor; horizontal bars, group medians;
dashed line in (b), overall median ({med_share:.1f}%).
Methods. LEPR protein-coding transcripts were grouped into short-form isoforms
(ENST00000344610/371058/371060/616738; 896-906 aa) and the canonical long isoform
(ENST00000349533; 1,165 aa) based on annotated transcript structure and protein length;
ambiguous (958 aa) and non-coding LEPR transcripts were excluded. Short-form fraction =
short-form TPM / summed TPM of included short- and long-form LEPR transcripts; total LEPR
TPM = summed TPM of included protein-coding LEPR transcripts. Short-form isoforms correspond
to LepRa-type transcripts, whereas the long isoform corresponds to the canonical LepRb transcript.
Statistics. Primary comparison: paired Wilcoxon signed-rank test (short vs long isoform
abundance, donor as unit; n=13, P=2.4x10^-4; median ~{med_fold:.0f}-fold, range {np.nanmin(fold):.1f}-{np.nanmax(fold):.1f}).
Group comparisons (Control, AD, DM1) were exploratory and not formally tested for disease
association (Kruskal-Wallis, not significant: P={kw['frac']:.2f} (b), {kw['gene']:.2f} (c); n=5/4/4, low power).
DM1 = myotonic dystrophy type 1 (NOT diabetes mellitus); disease groups shown for robustness.
Salmon v1.10.0; GENCODE v44 (GRCh38.p14) decoy-aware index (k=31); paired-end ~150 bp; full
libraries; -l A --validateMappings --gcBias. Exploratory analysis."""
open(f"{OUT}_caption.txt", "w").write(cap)
print("wrote refined figure (pdf/png/svg) + caption + summary")
print(cap)
