#!/usr/bin/env python3
"""
STEP 4 (OPTIONAL, NOT RUN) — cell-type deconvolution to estimate ChP epithelial
fraction per bulk sample. SKIPPED: requires a normal-human-ChP snRNA-seq reference
we do not currently hold (candidate: EMBO J 2024 normal human ChP; GEO accession TBD).

Core reviewer defense does NOT depend on this — STEP 3 (junction-vs-neuronal, null in
both cohorts) + STEP 1/2 (high TTR, no contamination-flagged sample) already stands.
Deconvolution is corroborative, not a blocker.

To activate when a reference is obtained:
  1. Drop reference into  ref/chp_snrna_reference.h5ad  (cells x genes, with cell_type).
  2. Pick a method (per-cohort; NEVER pool cohorts):
       - BisqueRNA / MuSiC (R, marker-free, robust to platform diff), OR
       - CIBERSORTx (signature matrix from the reference).
  3. Build a signature/reference profile from the snRNA-seq, restricted to genes
     shared with the bulk Salmon gene matrix (GENCODE v44 symbols).
  4. Deconvolve each cohort's bulk gene-TPM matrix -> per-sample fractions
     {epithelial, neuron, endothelial, astrocyte, immune, ...}.
  5. Deliverable: per-sample epithelial-fraction barplot + scatter of
     (neuron fraction) vs (LepRb junction CPM) — should mirror STEP 3 (no positive dep).

Inputs already available for this step:
  results/marker_matrix_GSE228458.csv , results/marker_matrix_GSE137619.csv
  (gene-level TPM; extend extract step to full-gene matrix if whole-transcriptome
   deconvolution input is needed — quant.sf on O2 scratch + tx2gene.tsv.)

Status: placeholder. No computation performed.
"""
import sys
print(__doc__)
print("STEP 4 deconvolution: SKIPPED (no reference). This is a placeholder; nothing was run.")
sys.exit(0)
