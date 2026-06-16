#!/usr/bin/env python3
"""Reproduce the leptin-receptor isoform-dominance result using the
`isoform-dominance` package as the analysis engine.

This script regenerates the key in-silico result of the study (short LEPR isoform
predominates over the long isoform in control human choroid plexus, across two
independent public cohorts) by calling the published `isoform-dominance` package
rather than any bespoke one-off code. It writes per-donor tables and a figure and
prints the paired Wilcoxon result.

Usage:
    pip install -r requirements.txt
    python reproduce.py

Expected: combined n = 11, 11/11 short > long, paired Wilcoxon P = 9.8e-4.
"""
import os
import csv
import tempfile

from isoform_dominance import extract, stats

# Per-transcript control choroid-plexus TPM (Salmon, GENCODE v44; GRCh38.p14).
# short LepRa = ENST00000371060 + ENST00000616738 ; long LepRb = ENST00000349533
DATA = {
    "GSE228458": [
        ("ctrl1", 6.521574, 36.556954, 1.314808), ("ctrl2", 1.250224, 8.637422, 0.349876),
        ("ctrl3", 2.238854, 13.496253, 0.704219), ("ctrl4", 3.684999, 15.015936, 0.721059),
        ("ctrl5", 0.495420, 6.599920, 1.440164),
    ],
    "GSE137619": [
        ("ctrl1", 9.130898, 9.898263, 0.351829), ("ctrl2", 5.896581, 5.928743, 0.216838),
        ("ctrl3", 7.516654, 20.317774, 1.033148), ("ctrl4", 1.804960, 3.059827, 0.143530),
        ("ctrl5", 3.730544, 2.848134, 0.203918), ("ctrl6", 2.582321, 4.337355, 0.166890),
    ],
}
CONFIG = {
    "gene": "LEPR",
    "groups": {"short": ["ENST00000371060", "ENST00000616738"], "long": ["ENST00000349533"]},
    "primary_comparison": ["short", "long"],
}


def _write_quant(path, t1, t2, tl):
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, "quant.sf"), "w") as f:
        f.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        f.write("ENST00000371060.5\t3000\t2800\t%s\t100\n" % t1)
        f.write("ENST00000616738.1\t3000\t2800\t%s\t100\n" % t2)
        f.write("ENST00000349533.11\t4000\t3800\t%s\t10\n" % tl)


def main():
    outdir = "results"
    os.makedirs(outdir, exist_ok=True)
    with tempfile.TemporaryDirectory() as base:
        perdonor = {}
        for cohort, rows in DATA.items():
            qd = os.path.join(base, cohort)
            for donor, t1, t2, tl in rows:
                _write_quant(os.path.join(qd, donor), t1, t2, tl)
            sm = os.path.join(base, "sample_map_%s.csv" % cohort)
            with open(sm, "w", newline="") as f:
                w = csv.writer(f); w.writerow(["donor", "condition"])
                for donor, *_ in rows:
                    w.writerow([donor, "control"])
            out = os.path.join(outdir, "perdonor_%s.csv" % cohort)
            extract.run(CONFIG, qd, sm, cohort, out)   # <- isoform-dominance does the TPM summation
            perdonor[cohort] = out

        res = stats.run(CONFIG, "control", perdonor, os.path.join(outdir, "lepr_dominance"))

    print("\nReproduced with isoform-dominance:")
    for name, n, ngt, p, fold in res["per_cohort"]:
        print("  %-12s n=%d  %d/%d short>long  fold=%.1fx  P=%.4g" % (name, n, ngt, n, fold, p))
    cn, cgt, cp, cfold = res["combined"]
    print("  COMBINED     n=%d  %d/%d  fold=%.1fx  P=%.4g" % (cn, cgt, cn, cfold, cp))
    print("\nFigure + tables written to ./results/")


if __name__ == "__main__":
    main()
