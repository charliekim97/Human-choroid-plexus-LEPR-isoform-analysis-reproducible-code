#!/usr/bin/env python3
"""
count_junctions.py — LEPR (and INSR control) isoform-usage junction counting.

Reproduces the splice-junction read-counting used for human ChP RNA-seq
(GSE228458). No aligner: equal-depth random subsample (seqkit, seed 42) ->
extract read sequences -> match isoform-defining junction k-mers with grep -F.

Method validated against Ensembl cDNA (each probe present only in its cognate
isoform; see scripts/junction_method_REPRODUCIBILITY.md).

LEPR last shared exon donor end = chr1:65,622,981 (GRCh38).
  LepRa (short, LEPR-205/208, acceptor 65,633,152)
  LepRb (long,  LEPR-202,     acceptor 65,636,191)

Usage:
  python3 count_junctions.py --fastq SRR2400xxxx_2.fastq.gz --target 30000000 \
      [--total <reads_in_file>]   # if --total omitted, p=1.0 (no subsample)

Requires: seqkit (>=2.x) on PATH; gzip/zcat. Writes counts to stdout.
"""
import argparse, subprocess, shlex, sys

SEED = 42

# --- isoform-defining junction probes (forward + reverse-complement) ---
def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

LEPRA = ["ATTTTCAGAAGAGAACGGAC"]            # short (LepRa); LEPR-205/208
LEPRB = ["ATTTTCAGAAGCCAGAAACG"]            # long  (LepRb); LEPR-202
INSR_INCL = ["TCGTCCCCAGAAAAACCTCT"]        # INSR exon10->11 (INSR-B), positive control
INSR_SKIP = ["TCGTCCCCAGGCCATCTCGG"]        # INSR exon10->12 (INSR-A skip)

def both(strands):
    out = []
    for s in strands:
        out += [s, rc(s)]
    return out

def write_patterns(path, kmers):
    with open(path, "w") as f:
        f.write("\n".join(kmers) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True, help="R2 FASTQ(.gz)")
    ap.add_argument("--target", type=int, default=30_000_000,
                    help="target reads for equal-depth subsample (default 30M)")
    ap.add_argument("--total", type=int, default=None,
                    help="total reads in file (to set subsample proportion p=target/total)")
    ap.add_argument("--mode", choices=["LEPR", "INSR"], default="LEPR")
    args = ap.parse_args()

    p = 1.0 if not args.total else min(1.0, args.target / args.total)

    if args.mode == "LEPR":
        short_pat, long_pat = both(LEPRA), both(LEPRB)
        slabel, llabel = "LepRa_short", "LepRb_long"
    else:
        short_pat, long_pat = both(INSR_INCL), both(INSR_SKIP)
        slabel, llabel = "INSR_include", "INSR_skip"
    write_patterns("/tmp/_short.txt", short_pat)
    write_patterns("/tmp/_long.txt", long_pat)
    write_patterns("/tmp/_all.txt", short_pat + long_pat)

    # seqkit sample (random, fixed seed) -> sequence lines -> reads overlapping either junction
    # single pass, no intermediate FASTQ written.
    cmd = (f"seqkit sample -s {SEED} -p {p:.6f} {shlex.quote(args.fastq)} "
           f"| awk 'NR%4==2' | grep -F -f /tmp/_all.txt")
    hits = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
    n_short = sum(1 for ln in hits.splitlines()
                  if any(k in ln for k in short_pat))
    n_long = sum(1 for ln in hits.splitlines()
                 if any(k in ln for k in long_pat))

    print(f"# fastq={args.fastq} mode={args.mode} seed={SEED} p={p:.4f} target={args.target}")
    print(f"{slabel}\t{n_short}")
    print(f"{llabel}\t{n_long}")
    if args.mode == "LEPR":
        tot = n_short + n_long
        frac = (n_long / tot) if tot else float("nan")
        print(f"LepRb_fraction\t{frac:.4f}")

if __name__ == "__main__":
    main()
