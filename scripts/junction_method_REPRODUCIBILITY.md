# LEPR isoform (LepRa vs LepRb) — human ChP bulk RNA-seq, reproducibility log
Dataset: GSE228458 / PRJNA949954 (Nutter 2023 Brain). 13 donors (5 control, 4 AD, 4 DM1). R2 reads.

## Diagnostic probes (junction-spanning, GRCh38)
LEPR last shared exon donor end chr1:65,622,981.
  LepRa (short, LEPR-205/208 accept 65,633,152):  ATTTTCAGAAGAGAACGGAC  (+rc GTCCGTTCTCTTCTGAAAAT)
  LepRb (long,  LEPR-202 accept 65,636,191):       ATTTTCAGAAGCCAGAAACG  (+rc CGTTTCTGGCTTCTGAAAAT)
Validated vs Ensembl cDNA: each junction present only in its own isoform (4/4 clean).

## Equal-depth random downsample (FINAL analysis)
sampling_method = seqkit sample (random), seed = 42, target ≈ 30,000,000 reads
per file proportion p = 30e6 / total_reads (single-pass streaming, no intermediate FASTQ)

Command pattern (per sample):
  seqkit sample -s 42 -p <p> <SRR>_2.fastq[.gz] \
    | awk 'NR%4==2' | grep -F -f JALL.txt > hits
  short = grep -F -f JS.txt hits | wc -l   # LepRa
  long  = grep -F -f JL.txt hits | wc -l   # LepRb
(JALL.txt = JS.txt + JL.txt. seqkit v2.8.2 linux_arm64.)

Runtime in-environment: gz 14-17s, raw 20GB ~28s (all < 43s timeout). No external server needed.

## Results
Final per-donor table: LEPR_isoform_FINAL.csv
  columns: donor,group,SRR,total_reads,sampling_method,seed,prop_sampled,n_reads_sampled,
           LepRa_reads,LepRb_reads,LepRa_CPM,LepRb_CPM,LepRb_fraction

Pooled (30M-normalized): control LepRa 245 : LepRb 6 | AD 171:6 | DM1 109:0
LepRb fraction: control 2.4%, AD 3.4%, DM1 0% -> LepRa dominant (>95%) in ALL groups.
control vs DM1 LepRa_CPM: Mann-Whitney p=0.46 (n=5 vs 4; DM1 trend lower but n.s.).

## Secondary (not yet run): INSR exon-11 positive control (DM1 spliceopathy)
  INSR_INCL TCGTCCCCAGAAAAACCTCT (+rc) ; INSR_SKIP TCGTCCCCAGGCCATCTCGG (+rc)
  validated vs INSR-B cDNA. Run same grep pattern if INSR check desired.

## Caveats
bulk tissue (not single-cell); short-read junction counts (not long-read);
n small (5/4/4); LepRb low-count -> fractions noisy. Direction (LepRa>>LepRb) robust.

## Prefix sanity check (NOT final; file-order biased) -> prefix_sanity.txt
