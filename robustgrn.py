#!/usr/bin/env python3
"""
RobustGRN: Mutual rank from CPM + Gene Regulatory Network inference via PyRRA.

Expects:
  - data/{genome_id}/cpm/*.cpm.tsv.gz   (input CPM)
  - tfs/{genome_id}.ptfdb.tsv           (TF list, first column = gene ID)

Produces:
  - data/{genome_id}/MR/{accession}.mr.tsv.gz  (mutual rank, from mutclust)
  - data/{genome_id}/{genome_id}.grn.tsv       (tf, target, p, adj; per-TF FDR, adj < 0.05 by default)

Usage:
  python robustgrn.py [--data-dir data] [--tfs-dir tfs] [--genome GENOME] [--mr-only | --grn-only]
"""

import os
import glob
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import numpy as np
from pyrra import aggregate_ranks
from statsmodels.stats.multitest import multipletests


# -----------------------------------------------------------------------------
# Defaults (mutclust mr: -m 200 -t 70 --log2; GRN: MR < 100, adj < 0.05)
# -----------------------------------------------------------------------------
MR_MUTCLUST_TOP = 200
MR_MUTCLUST_THREADS = 70
MR_MAX = 100
ADJ_CUTOFF_DEFAULT = 0.05
TSV_SEP = "\t"
MIN_SAMPLES = 12


def _cpm_sample_count(cpm_path):
    """Return number of columns (samples) in first data row; first column is gene ID."""
    if cpm_path.endswith(".gz"):
        import gzip
        with gzip.open(cpm_path, "rt") as f:
            line = f.readline()
    else:
        with open(cpm_path) as f:
            line = f.readline()
    return len(line.strip().split(TSV_SEP)) - 1


def run_mr_for_cpm(cpm_path, mr_path, top=MR_MUTCLUST_TOP, threads=MR_MUTCLUST_THREADS, skip_existing=True):
    """Run mutclust mr on one CPM file; write MR to mr_path. Returns True if ran or skipped (exists), False if skipped (low samples)."""
    if skip_existing and os.path.exists(mr_path):
        return True
    n_cols = _cpm_sample_count(cpm_path)
    if n_cols < MIN_SAMPLES:
        return False
    os.makedirs(os.path.dirname(mr_path), exist_ok=True)
    subprocess.run(
        ["mutclust", "mr", "--input", cpm_path, "--output", mr_path, "-m", str(top), "-t", str(threads), "--log2"],
        check=True,
    )
    return True


def run_mr_for_genome(genome_id, data_dir="data", top=MR_MUTCLUST_TOP, threads=MR_MUTCLUST_THREADS, skip_existing=True):
    """For one genome, run mutclust mr on each CPM in data/{genome_id}/cpm/; write to data/{genome_id}/MR/{accession}.mr.tsv.gz."""
    cpm_dir = os.path.join(data_dir, genome_id, "cpm")
    mr_dir = os.path.join(data_dir, genome_id, "MR")
    if not os.path.isdir(cpm_dir):
        return 0, 0
    cpm_files = sorted(glob.glob(os.path.join(cpm_dir, "*.cpm.tsv.gz")))
    ran, skipped = 0, 0
    for cpm_path in cpm_files:
        accession = os.path.basename(cpm_path).replace(".cpm.tsv.gz", "")
        mr_path = os.path.join(mr_dir, f"{accession}.mr.tsv.gz")
        ok = run_mr_for_cpm(cpm_path, mr_path, top=top, threads=threads, skip_existing=skip_existing)
        if ok:
            if os.path.exists(mr_path) and not skip_existing:
                ran += 1
            else:
                skipped += 1
        else:
            skipped += 1
    return ran, skipped


def load_tfs(genome_id, tfs_dir="tfs"):
    """Load TF gene IDs from tfs/{genome_id}.ptfdb.tsv (first column)."""
    path = os.path.join(tfs_dir, f"{genome_id}.ptfdb.tsv")
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path, sep=TSV_SEP, header=None, usecols=[0])
    return df[0].astype(str).tolist()


def get_mr_files(genome_id, data_dir="data"):
    """Return sorted list of MR files for this genome: data/{genome_id}/MR/*.mr.tsv.gz."""
    mr_dir = os.path.join(data_dir, genome_id, "MR")
    pattern = os.path.join(mr_dir, "*.mr.tsv.gz")
    return sorted(glob.glob(pattern))


def parse_mr_file_for_tfs(mr_path, tfs_set, mr_max=MR_MAX):
    """Read one MR file; filter MR < mr_max; return dict tf -> list of (partner, mr) for this file. Uses pigz -dc if .gz."""
    if mr_path.endswith(".gz"):
        try:
            pigz_proc = subprocess.Popen(["pigz", "-dc", mr_path], stdout=subprocess.PIPE)
            df = pd.read_csv(pigz_proc.stdout, sep=TSV_SEP, usecols=["Gene1", "Gene2", "MR"])
            pigz_proc.wait()
        except FileNotFoundError:
            import gzip
            df = pd.read_csv(mr_path, sep=TSV_SEP, usecols=["Gene1", "Gene2", "MR"], compression="gzip")
    else:
        df = pd.read_csv(mr_path, sep=TSV_SEP, usecols=["Gene1", "Gene2", "MR"])
    df = df[df["MR"] < mr_max]
    mask1 = df["Gene1"].isin(tfs_set) & (df["Gene2"] != df["Gene1"])
    mask2 = df["Gene2"].isin(tfs_set) & (df["Gene1"] != df["Gene2"])
    pairs1 = df.loc[mask1, ["Gene1", "Gene2", "MR"]].rename(columns={"Gene1": "tf", "Gene2": "partner"})
    pairs2 = df.loc[mask2, ["Gene1", "Gene2", "MR"]].rename(columns={"Gene2": "tf", "Gene1": "partner"})
    combined = pd.concat([pairs1[["tf", "partner", "MR"]], pairs2[["tf", "partner", "MR"]]], ignore_index=True)
    combined = combined.sort_values(["tf", "MR"])
    tf_partners = {
        tf: list(zip(grp["partner"].tolist(), grp["MR"].tolist()))
        for tf, grp in combined.groupby("tf", sort=False)
    }
    return tf_partners


def ranked_list_for_tf(partner_mr_list):
    """Turn list of (partner, mr) into list of gene IDs (already sorted by MR)."""
    if not partner_mr_list:
        return []
    return [g for g, _ in partner_mr_list]


def build_per_tf_ranked_lists(genome_id, tfs_list, data_dir="data", mr_max=MR_MAX):
    """For each MR file, for each TF: append that TF's ranked list from this file. Returns (tf -> list of lists), n_mr_files."""
    tfs_set = set(tfs_list)
    mr_files = get_mr_files(genome_id, data_dir=data_dir)
    if not mr_files:
        return {}, 0
    per_tf_ranked_lists = {tf: [] for tf in tfs_list}
    for mr_path in mr_files:
        tf_partners = parse_mr_file_for_tfs(mr_path, tfs_set, mr_max=mr_max)
        for tf in tfs_list:
            ranked = ranked_list_for_tf(tf_partners.get(tf, []))
            per_tf_ranked_lists[tf].append(ranked)
    return per_tf_ranked_lists, len(mr_files)


def _pyrra_one_tf(args):
    tf, glists, method = args
    if len(glists) < 1:
        return tf, []
    try:
        agg = aggregate_ranks(glist=glists, method=method)
    except Exception:
        return tf, []
    if agg is None or agg.empty:
        return tf, []
    return tf, [{"tf": tf, "target": r["Name"], "p": r["Score"]} for _, r in agg.iterrows()]


def run_pyrra_per_tf(per_tf_ranked_lists, method="RRA", n_jobs=1):
    """Run PyRRA on each TF's list of ranked lists; return DataFrame (tf, target, p)."""
    items = []
    for tf, glists in per_tf_ranked_lists.items():
        non_empty = [g for g in glists if g]
        if len(non_empty) >= 1:
            items.append((tf, non_empty, method))
    if not items:
        return pd.DataFrame(columns=["tf", "target", "p", "adj"])
    if n_jobs <= 1:
        rows = []
        for item in items:
            _, r = _pyrra_one_tf(item)
            rows.extend(r)
    else:
        rows = []
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            for future in as_completed(executor.submit(_pyrra_one_tf, item) for item in items):
                _, r = future.result()
                rows.extend(r)
    if not rows:
        return pd.DataFrame(columns=["tf", "target", "p", "adj"])
    return pd.DataFrame(rows)


def fdr_per_tf(df):
    """Per-TF FDR (Benjamini–Hochberg). Adds column 'adj'."""
    adj = np.full(len(df), np.nan)
    for tf, idx in df.groupby("tf", sort=False).indices.items():
        p = df.loc[idx, "p"].values
        _, p_adj, _, _ = multipletests(p, method="fdr_bh")
        adj[idx] = p_adj
    df = df.copy()
    df["adj"] = adj
    return df


def run_grn_for_genome(genome_id, data_dir="data", tfs_dir="tfs", mr_max=MR_MAX, method="RRA", adj_cutoff=ADJ_CUTOFF_DEFAULT, skip_existing=True, n_jobs=1):
    """Build GRN for one genome: load TFs, build ranked lists from data/{g}/MR/, PyRRA, per-TF FDR, write data/{g}/{g}.grn.tsv."""
    out_path = os.path.join(data_dir, genome_id, f"{genome_id}.grn.tsv")
    if skip_existing and os.path.exists(out_path):
        print(f"  -> Skip (exists): {out_path}")
        return
    tfs_list = load_tfs(genome_id, tfs_dir=tfs_dir)
    if not tfs_list:
        print(f"  -> No TFs for {genome_id}")
        return
    per_tf_ranked_lists, n_mr_files = build_per_tf_ranked_lists(genome_id, tfs_list, data_dir=data_dir, mr_max=mr_max)
    n_tfs_with_data = sum(1 for glists in per_tf_ranked_lists.values() if any(g for g in glists))
    if n_tfs_with_data == 0:
        print(f"  -> No MR data for any TF in {genome_id}")
        return
    print(f"  -> Loaded {n_mr_files} MR files; running PyRRA on {n_tfs_with_data} TFs...")
    df = run_pyrra_per_tf(per_tf_ranked_lists, method=method, n_jobs=n_jobs)
    if df.empty:
        print(f"  -> No PyRRA results for {genome_id}")
        return
    df = fdr_per_tf(df)
    df = df[df["adj"] < adj_cutoff]
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.to_csv(out_path, sep=TSV_SEP, index=False)
    print(f"  -> {out_path}: {len(df)} edges (adj < {adj_cutoff}), {n_tfs_with_data} TFs, {n_mr_files} MR files")


def main():
    parser = argparse.ArgumentParser(description="RobustGRN: MR from CPM + GRN via PyRRA (per-TF FDR).")
    parser.add_argument("--data-dir", default="data", help="Data root: data/{genome_id}/cpm/, data/{genome_id}/MR/, data/{genome_id}/{genome_id}.grn.tsv")
    parser.add_argument("--tfs-dir", default="tfs", help="TF lists: tfs/{genome_id}.ptfdb.tsv")
    parser.add_argument("--genome", default=None, help="Process only this genome_id")
    parser.add_argument("--mr-only", action="store_true", help="Only run mutclust MR (skip GRN)")
    parser.add_argument("--grn-only", action="store_true", help="Only run GRN (skip MR; assume MR files exist)")
    parser.add_argument("--mr-max", type=float, default=MR_MAX, help=f"Keep gene pairs with MR < this for GRN (default: {MR_MAX})")
    parser.add_argument("--method", default="RRA", choices=["RRA", "stuart", "min", "mean", "median", "geom.mean"], help="PyRRA method")
    parser.add_argument("--threads", type=int, default=None, metavar="N", help="PyRRA workers per genome (0 = all cores)")
    parser.add_argument("--adj-cutoff", type=float, default=ADJ_CUTOFF_DEFAULT, metavar="F", help=f"Keep only edges with FDR-adjusted p-value (adj) < this (default: {ADJ_CUTOFF_DEFAULT})")
    parser.add_argument("--no-skip-existing", action="store_true", help="Recompute existing MR/GRN outputs")
    args = parser.parse_args()

    skip_existing = not args.no_skip_existing
    if args.threads is None or args.threads == 1:
        n_jobs = 1
    elif args.threads == 0:
        n_jobs = max(1, os.cpu_count() or 1)
    else:
        n_jobs = max(1, args.threads)

    # Resolve genomes from data dir
    if args.genome:
        genome_ids = [args.genome]
    else:
        genome_ids = []
        if os.path.isdir(args.data_dir):
            for name in os.listdir(args.data_dir):
                path = os.path.join(args.data_dir, name)
                if os.path.isdir(path) and os.path.isdir(os.path.join(path, "cpm")):
                    genome_ids.append(name)
        genome_ids = sorted(genome_ids)

    if not genome_ids:
        print("No genomes found (expect data/{genome_id}/cpm/).")
        return

    # Step 1: MR from CPM (unless --grn-only)
    if not args.grn_only:
        print("=== Step 1: Mutual rank (mutclust mr) ===")
        for genome_id in genome_ids:
            cpm_dir = os.path.join(args.data_dir, genome_id, "cpm")
            if not os.path.isdir(cpm_dir):
                continue
            print(f"--- {genome_id} ---")
            for cpm_path in sorted(glob.glob(os.path.join(cpm_dir, "*.cpm.tsv.gz"))):
                accession = os.path.basename(cpm_path).replace(".cpm.tsv.gz", "")
                mr_path = os.path.join(args.data_dir, genome_id, "MR", f"{accession}.mr.tsv.gz")
                if skip_existing and os.path.exists(mr_path):
                    print(f"  Skip: {mr_path}")
                    continue
                n_cols = _cpm_sample_count(cpm_path)
                if n_cols < MIN_SAMPLES:
                    print(f"  Skip (samples < {MIN_SAMPLES}): {cpm_path}")
                    continue
                print(f"  MR: {cpm_path} -> {mr_path}")
                os.makedirs(os.path.dirname(mr_path), exist_ok=True)
                subprocess.run(
                    ["mutclust", "mr", "--input", cpm_path, "--output", mr_path, "-m", str(MR_MUTCLUST_TOP), "-t", str(MR_MUTCLUST_THREADS), "--log2"],
                    check=True,
                )

    # Step 2: GRN (unless --mr-only)
    if not args.mr_only:
        grn_genomes = [g for g in genome_ids if os.path.isdir(os.path.join(args.data_dir, g, "MR")) and os.path.exists(os.path.join(args.tfs_dir, f"{g}.ptfdb.tsv"))]
        if not grn_genomes:
            print("No genomes with both data/{g}/MR/ and tfs/{g}.ptfdb.tsv for GRN.")
        else:
            print(f"=== Step 2: GRN (PyRRA + per-TF FDR), adj < {args.adj_cutoff} ===")
            for genome_id in grn_genomes:
                print(f"--- {genome_id} ---")
                run_grn_for_genome(
                    genome_id,
                    data_dir=args.data_dir,
                    tfs_dir=args.tfs_dir,
                    mr_max=args.mr_max,
                    method=args.method,
                    adj_cutoff=args.adj_cutoff,
                    skip_existing=skip_existing,
                    n_jobs=n_jobs,
                )


if __name__ == "__main__":
    main()
