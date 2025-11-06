# GRNclust

A streamlined and reproducible Snakemake workflow for batch-correcting RNA-seq gene expression data (Combat-Seq), generating Gene Regulatory Networks (GRNs; GRNBoost2), and clustering GRNs with the Leiden algorithm.

---

## Installation

### Using Conda
1. Clone the repository:
   ```
   bash
   git clone https://github.com/eporetsky/GRNclust.git
   cd GRNclust
   ```

2. Create the Conda environment: 
    ```
    conda env create -f environment.yml
    conda activate grnclust
    
    # fixes an issue with GRNboost2
    bash post-install.sh 
    
    # arboreto core.py uses 1 thread by default, replace with N
    bash post-threads.sh N
    ```

* Note: If you still get `TypeError: Must supply at least one delayed object` with GRNboost2, read the my comments in `post-install.sh`.

---

## Workflow Overview

The workflow consists of the following steps:

1. **Preprocessing**:
   - Merges multiple bulk RNA-seq count files into a single dataframe.
   - Combines and batch-corrects gene expression data from multiple input files.
   - Outputs a normalized and log-scaled expression matrix.

2. **GRN Inference**:
   - Uses the `GRNBoost2` algorithm to infer a gene regulatory network (GRN) from the processed expression data and a list of transcription factors (TFs).

3. **Leiden Clustering**:
   - Applies the Leiden clustering algorithm to the inferred GRN.
   - Filters edges to retain only those within the same cluster.

---

## Command to Run the workflow

To execute the workflow on multiple count files with batch correction use the following command:

```bash
python -m snakemake --cores 1 --config input_folder=example output_folder=example tf_file=example/AtCol-0.ptfdb.tsv resolution=0.1 threads=12
```

## Run GRN analysis on a single expression file

```bash
python scripts/run_grnboost2.py --expression_file acc.cpm.tsv.gz \
                                --tf_file genotype.ptfdb.tsv \
                                --output_file acc.grn.tsv
````

---

## Workflow Details

### **1. Preprocessing**

- **Input**:
  - Dynamically detects all count (tsv) files in the `counts` subdirectory of the `input_folder`.
  - Example: `example/counts/sample1.counts.tsv.gz`, `example/counts/sample2.counts.tsv.gz`.

- **Output**:
  - A batch-corrected and normalized expression matrix: `example/expression.batch.logscale.tsv`.

- **Script**:
  - The preprocessing is performed using the `preprocess.R` script, which:
    - Merges input files by `geneID`.
    - Performs batch correction using `ComBat-Seq`.
    - Log-transforms and scales the data.

---

### **2. GRN Inference**

- **Input**:
  - The batch-corrected expression matrix: `example/expression.batch.logscale.tsv`.
  - A list of transcription factors: `example/AtCol-0.ptfdb.tsv`.

- **Output**:
  - The inferred GRN: `example/complete_grn.csv`.

- **Script**:
  - The GRN inference is performed using the `GRNBoost2` algorithm in the `run_grnboost2.py` script.

---

### **3. Leiden Clustering**

- **Input**:
  - The inferred GRN: `example/complete_grn.csv`.

- **Output**:
  - A filtered GRN with edges retained only within clusters: `example/filtered_grn_within_clusters_res0.1.csv`.

- **Script**:
  - The clustering is performed using the Leiden algorithm in the `run_leiden_clustering.py` script.

---

## Command-Line Arguments

The workflow accepts the following arguments via the `--config` option:

| Argument         | Description                                                                 |
|-------------------|-----------------------------------------------------------------------------|
| `input_folder`    | Path to the folder containing the input files (e.g., `example`).           |
| `output_folder`   | Path to the folder where output files will be saved (e.g., `example`).     |
| `tf_file`         | Path to the transcription factor list file (e.g., `example/AtCol-0.ptfdb.tsv`). |
| `resolution`      | Resolution parameter for Leiden clustering (e.g., `0.1`).                 |
| `threads`         | Number of threads to use for parallel processing (e.g., `12`).            |

---

### Notes

- The workflow dynamically detects input files in the `counts` subdirectory of the `input_folder`.
- The `resolution` parameter allows you to adjust the granularity of the Leiden clustering.
- The workflow tries to avoid unnecessary re-runs by properly tracking input-output dependencies.

---

## Citation

This workflow uses the following software. Please consider including the relecant citations:
- [ComBat-seq](https://github.com/zhangyuqing/ComBat-seq)
- [GRNboost2](https://github.com/aertslab/GRNBoost) - available via the [3-Clause BSD license](https://opensource.org/licenses/BSD-3-Clause).
- [leidenalg](https://leidenalg.readthedocs.io/en/stable/index.html)
- [PlantTFDB](https://planttfdb.gao-lab.org/) - used for predicting TFs in the example file