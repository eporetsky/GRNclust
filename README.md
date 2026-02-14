# RobustGRN

**RobustGRN** builds gene regulatory networks (GRNs) from RNA-seq expression data by (1) computing Mutual Rank (MR) co-expression from CPM tables using [MutClust](https://github.com/eporetsky/mutclust), and (2) aggregating data for individual transcription factors over multiple experiments using [PyRRA](https://github.com/eporetsky/pyrra) (Robust Rank Aggregation). It is meant as a **meta-analysis GRN inference** based on **Mutual-Rank coexpression networks**.

The method was developed for **PlantApp**; results for multiple plant species are available at [PlantApp GRN](https://www.plantapp.org/apps/grn). 

### Disclaimer

This method is experimental and has not been benchmarked against existing GRN inference methods nor validated against experimental data (e.g. ChIP-seq, perturbation assays). There is no current plan for any further active development. The method prioritizes speed and efficiency with a goal of identifying strong connections between TFs and target genes. Results are co-expression–based and do not imply direct regulation. Users are responsible for evaluating and interpreting outputs for their own applications. For more established GRN methods, see [GRNbenchmark public methods](https://grnbenchmark.org/public-methods) or other resources.

## Data layout

All inputs and outputs use a single root directory (default `data/`):

```
data/
  {genome_id}/
    cpm/                          # Input: CPM tables (genes × samples)
      {accession}.cpm.tsv.gz
    MR/                           # Output: mutual rank (from mutclust)
      {accession}.mr.tsv.gz
    {genome_id}.grn.tsv           # Output: GRN edges (tf, target, p, adj)
tfs/
  {genome_id}.ptfdb.tsv           # Input: TF list (first column = gene ID)
```

- **CPM**: one gzipped TSV per experiment; first column = gene ID, remaining columns = samples (CPM).
- **MR**: gene–gene mutual rank tables produced by `mutclust mr` (Gene1, Gene2, MR).
- **TF list**: tab-separated, no header; column 1 = TF gene ID, column 2 = family (optional). [PlantTFDB](https://planttfdb.gao-lab.org) can be used for predicting TFs.

## Requirements

- **Python 3** with: `pandas`, `numpy`, `pyrra`, `statsmodels`
- **mutclust** (CLI) for the MR step: [mutclust](https://github.com/eporetsky/mutclust)
- **pigz** (optional): speeds up reading gzipped MR files

### Install Python dependencies

**Conda (from environment file):**

```bash
conda env create -f environment.yml
conda activate robustgrn
```

## Usage

Run both steps (MR from CPM, then GRN) for all genomes under `data/`:

```bash
python robustgrn.py --data-dir data --tfs-dir tfs --adj-cutoff 1 --threads 72
```

Options:

| Option | Default | Description |
|--------|---------|-------------|
| `--data-dir` | `data` | Root for `{genome_id}/cpm/`, `{genome_id}/MR/`, `{genome_id}/{genome_id}.grn.tsv` |
| `--tfs-dir` | `tfs` | Directory of `{genome_id}.ptfdb.tsv` TF lists |
| `--genome` | (all) | Process only this genome ID |
| `--mr-only` | off | Only run mutclust MR; do not run GRN |
| `--grn-only` | off | Only run GRN; assume MR files already exist |
| `--mr-max` | 100 | Use only gene pairs with MR &lt; this in GRN |
| `--method` | RRA | PyRRA method: RRA, stuart, min, mean, median, geom.mean |
| `--threads` | 1 | Parallel workers for PyRRA per genome (0 = all cores) |
| `--adj-cutoff` | 0.05 | Keep only edges with FDR-adjusted p-value (adj) &lt; this |
| `--no-skip-existing` | off | Recompute existing MR/GRN outputs |

Examples:

```bash
# Single genome, 8 threads for PyRRA
python robustgrn.py --genome ZmB73 --threads 8
```

## Pipeline summary

1. **MR (mutclust)**  
   For each `data/{genome_id}/cpm/*.cpm.tsv.gz`: run `mutclust mr` with `-m 200 -t 70 --log2`. Experiments with &lt; 12 samples are skipped. Output: `data/{genome_id}/MR/{accession}.mr.tsv.gz`.

2. **GRN (PyRRA + FDR)**  
   - Load TFs from `tfs/{genome_id}.ptfdb.tsv`.
   - For each MR file in `data/{genome_id}/MR/*.mr.tsv.gz`, for each TF: extract partners with MR &lt; `--mr-max`, rank by MR → one ranked list per (TF, experiment).
   - Run **PyRRA** on each TF’s list of ranked lists → one p-value per (TF, target).
   - Apply **per-TF FDR** (Benjamini–Hochberg) for the `adj` column.
   - Keep only edges with **adj &lt; --adj-cutoff** (default 0.05); write `data/{genome_id}/{genome_id}.grn.tsv` with columns: `tf`, `target`, `p`, `adj`.

## Output: GRN TSV

`data/{genome_id}/{genome_id}.grn.tsv`:

- **tf** – transcription factor (gene ID)
- **target** – target gene ID
- **p** – raw p-value from PyRRA
- **adj** – FDR-adjusted p-value (per TF)

Only rows with **adj** below the `--adj-cutoff` threshold (default 0.05) are written.