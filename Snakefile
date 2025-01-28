import os
import glob

# Dynamically find all count files in the input folder
def get_count_files():
    return glob.glob(f"{config['input_folder']}/counts/*")

rule all:
    input:
        f"{config['output_folder']}/filtered_grn_within_clusters_res{config['resolution']}.csv"

rule preprocess:
    input:
        counts=get_count_files()
    output:
        f"{config['input_folder']}/expression.batch.logscale.tsv"
    threads: config["threads"]
    shell:
        "Rscript scripts/preprocess.R {config[input_folder]} {config[input_folder]}"

rule grn_inference:
    input:
        expression=f"{config['input_folder']}/expression.batch.logscale.tsv",
        tf_list=config["tf_file"]
    output:
        f"{config['output_folder']}/complete_grn.csv"
    threads: config["threads"]
    shell:
        "python scripts/run_grnboost2.py --expression_file {input.expression} --tf_file {input.tf_list} --output_file {output}"

rule leiden_clustering:
    input:
        grn=f"{config['output_folder']}/complete_grn.csv"
    output:
        f"{config['output_folder']}/filtered_grn_within_clusters_res{config['resolution']}.csv"
    params:
        resolution=config["resolution"]
    threads: config["threads"]
    shell:
        "python scripts/run_leiden_clustering.py --grn_file {input.grn} --resolution {params.resolution} --output_file {output}"