# Human choroid plexus LEPR isoform analysis

Reproducible code for the in-silico human-relevance analysis of LEPR isoform usage in
human choroid plexus (ChP), supporting the manuscript on a leptin receptor–LRP1 complex
in the choroid plexus and tanycytes.

Two independent public human ChP bulk RNA-seq datasets were re-analysed from raw reads with
two complementary read-outs — transcript-level isoform **expression** (Salmon) and
splice-junction isoform **usage** — plus a marker-based choroid-plexus contamination control.
The two cohorts were processed with the same workflow but analysed as separate discovery and
replication cohorts (never pooled).

## Datasets (public)

| Role | GEO | Study | Tissue | n |
|------|-----|-------|--------|---|
| Discovery | GSE228458 | Nutter et al., *Brain* 2023 | post-mortem human ChP, bulk | 5 control / 4 AD / 4 DM1 |
| Replication | GSE137619 | Rodríguez-Lorenzo et al., *Acta Neuropathol Commun* 2020 | post-mortem human ChP, bulk | 6 control / 6 MS |

FASTQ were obtained from the ENA. DM1 = myotonic dystrophy type 1 (not diabetes mellitus).
Controls are neurologically unaffected / non-neurological post-mortem donors.

## Software

| Tool | Version |
|------|---------|
| Salmon | 1.10.0 (precompiled linux x86_64) |
| Reference | GENCODE v44 (GRCh38.p14), decoy-aware index, k=31 |
| seqkit | 2.8.2 |
| Python | 3.10 (numpy, scipy, matplotlib, openpyxl) |

See `requirements.txt`.

## Repository structure

```
scripts/    analysis pipeline (see order below)
data/       per-donor result tables produced by the pipeline (LEPR TPM, junction counts, marker matrices, scores)
```

## Pipeline / order to reproduce

1. **Quantification (HPC).** `scripts/salmon_LEPR_O2.sbatch`, `scripts/salmon_GSE137619_O2.sbatch`
   — download FASTQ from ENA, build the decoy-aware GENCODE v44 index, run `salmon quant`
   (`-l A --validateMappings --gcBias`), per cohort.
2. **LEPR isoform expression.** `scripts/extract_lepr.py` — per-donor LEPR transcript TPM,
   short (LepRa-type) vs long (LepRb, ENST00000349533) grouping.
3. **Splice-junction usage.** `scripts/count_junctions.py` — equal-depth subsample (seqkit,
   seed 42) and exact junction-probe counting (see `scripts/junction_method_REPRODUCIBILITY.md`).
4. **Isoform-definition variants / figures.** `reanalyze_2tx.py` (2-transcript LepRa),
   `plot_individual_vs_long.py` (per-transcript), `plot_*` (display figures; final figures
   were redrawn in GraphPad Prism).
5. **Choroid-plexus contamination control.** `scripts/extract_markers.py` (marker-gene TPM
   from quant.sf) then `scripts/contamination_qc.py` (purity score, MAD outlier flag,
   LEPR-long vs neuronal-score Spearman correlation per cohort).
6. **(optional) deconvolution** — `scripts/step4_deconvolution_PLACEHOLDER.py` (not run; no reference).

## Outputs (in `data/`)

Per-donor LEPR transcript TPM, short/long summaries, junction-spanning counts (depth-normalized),
marker-gene TPM matrices, and per-sample contamination scores/flags, for both cohorts.

## Notes

- No new algorithms were developed; all third-party software is cited above with versions.
- Salmon quantification was run on the Harvard Medical School O2 high-performance cluster.
- Statistics: donor-level two-sided Wilcoxon signed-rank (`scipy.stats.wilcoxon`, exact);
  contamination control by Spearman correlation.

## Reproduce via isoform-dominance
pip install -r requirements.txt
python reproduce.py

## Citation

If you use this code, please cite the associated manuscript (details to be added on publication).

## License

MIT (see `LICENSE`).
