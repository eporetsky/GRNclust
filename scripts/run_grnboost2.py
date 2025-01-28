import pandas as pd
from arboreto.algo import grnboost2
import argparse

def run_grnboost2(expression_file, tf_file, output_file):
    # Load transcription factors
    tf_df = pd.read_csv(tf_file, sep="\t", header=None)
    tf_names = tf_df[0].tolist()

    # Load expression data
    data_corrected = pd.read_csv(expression_file, sep="\t")
    data_corrected = data_corrected.set_index("geneID")  # Set geneID as the index

    # Perform GRN inference
    grn_data = grnboost2(expression_data=data_corrected.T, tf_names=tf_names)

    # Save the GRN to a file
    grn_data.to_csv(output_file, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run GRNBoost2 to infer the gene regulatory network")
    parser.add_argument("--expression_file", required=True, help="Path to the batch-corrected and normalized expression file")
    parser.add_argument("--tf_file", required=True, help="Path to the transcription factor list file")
    parser.add_argument("--output_file", required=True, help="Path to save the inferred GRN")
    args = parser.parse_args()

    run_grnboost2(args.expression_file, args.tf_file, args.output_file)
