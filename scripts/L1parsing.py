#!/usr/bin/env python3
"""
rmsk_l1_div_plot.py
- Reads a standard RepeatMasker .out file (fixed-width layout with dashed header)
- Filters LINE-1 (L1) records
- Plots histogram of % divergence and saves CSV + summary

Usage:
  python rmsk_l1_div_plot.py --input hs1.repeatMasker.out --out hs1_L1 --bins 40 --top-families 0
"""

import argparse, re, csv, sys
from statistics import mean, median, pstdev
import matplotlib.pyplot as plt

L1_CLASS_PAT = re.compile(r"LINE\W*/*\W*L1", re.IGNORECASE)
L1_FAM_PAT = re.compile(r"(L1HS|L1PA\d+|L1PB\d+|L1[A-Z]{1,2}\d*|L1P\w+)", re.IGNORECASE)

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input","-i", required=True, help=".out file")
    ap.add_argument("--out","-o", default="l1_divergence", help="output prefix")
    ap.add_argument("--bins", type=int, default=40, help="histogram bins")
    ap.add_argument("--top-families", type=int, default=0, help="overlay top N families (0 = skip)")
    return ap.parse_args()

def is_l1(repClassFamily, repName):
    s = f"{repClassFamily} {repName or ''}"
    if L1_CLASS_PAT.search(s): return True
    if repName and repName.upper().startswith("L1"): return True
    return False

def infer_family(repName, repClassFamily):
    # Prefer explicit family after slash, else guess from repName, else Unknown
    fam = ""
    if "/" in (repClassFamily or ""):
        parts = repClassFamily.split("/", 1)
        fam = parts[1]
    if not fam and repName:
        m = L1_FAM_PAT.search(repName)
        if m: fam = m.group(0)
    return fam or "Unknown"

def read_rmsk_out(path):
    # find dashed line (column header underline)
    skip = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for i, line in enumerate(fh):
            if re.match(r"^\s*-+\s*$", line):
                skip = i + 1
                break
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for i, line in enumerate(fh):
            if i < skip: continue
            if not line.strip(): continue
            parts = re.split(r"\s+", line.strip())
            # Expected >= 14 columns from .out after header
            if len(parts) < 14: continue
            try:
                perc_div = float(parts[1])         # "% div."
            except ValueError:
                continue
            query   = parts[4]                     # sequence (chrom)
            q_begin = parts[5].replace(",", "")
            q_end   = parts[6].replace(",", "")
            try:
                start = int(q_begin); end = int(q_end)
            except ValueError:
                start = end = None
            strand = "-" if "C" in parts[8:10] else "+"
            repName = parts[9]                     # matching repeat name (e.g., L1HS)
            class_family = parts[10]               # e.g., LINE/L1

            if not is_l1(class_family, repName): 
                continue

            rows.append({
                "seqname": query,
                "start": start,
                "end": end,
                "strand": strand,
                "repName": repName,
                "repClassFamily": class_family,
                "repFamily": infer_family(repName, class_family),
                "perc_div": perc_div
            })
    return rows

def write_outputs(rows, out_prefix, bins, top_families):
    if not rows:
        print("No LINE-1 rows found with usable % divergence."); sys.exit(1)

    # CSV
    with open(f"{out_prefix}.csv","w",newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"Saved CSV → {out_prefix}.csv  (n={len(rows)})")

    vals = [r["perc_div"] for r in rows]
    n = len(vals); mn = mean(vals); md = median(vals)
    sd = pstdev(vals) if n > 1 else 0.0; mnv = min(vals); mxv = max(vals)
    with open(f"{out_prefix}_summary.txt","w") as f:
        f.write(f"n_L1: {n}\nmean_div: {mn:.4f}\nmedian_div: {md:.4f}\n"
                f"stdev_div: {sd:.4f}\nmin_div: {mnv:.4f}\nmax_div: {mxv:.4f}\n")
    print(f"Saved summary → {out_prefix}_summary.txt")

    # Main histogram
    plt.figure()
    plt.hist(vals, bins=bins)
    plt.xlabel("% divergence"); plt.ylabel("Count of LINE-1 hits")
    plt.title("RepeatMasker LINE-1 divergence")
    plt.tight_layout(); plt.savefig(f"{out_prefix}.png")
    print(f"Saved plot → {out_prefix}.png")

    # Optional: overlay by family
    if top_families and top_families > 0:
        from collections import Counter, defaultdict
        fam_counts = Counter([r["repFamily"] for r in rows])
        top = [fam for fam,_ in fam_counts.most_common(top_families)]
        plt.figure()
        for fam in top:
            sub = [r["perc_div"] for r in rows if r["repFamily"] == fam]
            if sub: plt.hist(sub, bins=bins, alpha=0.5, label=fam)
        plt.xlabel("% divergence"); plt.ylabel("Count")
        plt.title(f"Top {len(top)} L1 families: divergence")
        plt.legend(); plt.tight_layout()
        plt.savefig(f"{out_prefix}_families.png")
        print(f"Saved family overlay → {out_prefix}_families.png")

def main():
    args = parse_args()
    rows = read_rmsk_out(args.input)
    write_outputs(rows, args.out, args.bins, args.top_families)

if __name__ == "__main__":
    main()
