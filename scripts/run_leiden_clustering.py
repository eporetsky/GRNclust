import pandas as pd
import igraph as ig
import leidenalg as la
import argparse

def run_leiden_clustering(grn_file, resolution, output_file):
    # Load the GRN
    grn_data = pd.read_csv(grn_file)

    # Create an igraph graph from the GRN data
    edges = list(zip(grn_data['TF'], grn_data['target']))
    weights = grn_data['importance'].tolist()
    graph = ig.Graph.TupleList(edges, directed=True)
    graph.es['weight'] = weights

    # Apply Leiden clustering
    partition = la.find_partition(graph, la.CPMVertexPartition, resolution_parameter=resolution, weights=graph.es['weight'])

    # Extract cluster membership
    clusters = partition.membership
    graph.vs['cluster'] = clusters

    # Remove edges between nodes in different clusters
    edges_to_remove = [e.index for e in graph.es if graph.vs[e.source]['cluster'] != graph.vs[e.target]['cluster']]
    graph.delete_edges(edges_to_remove)

    # Save the filtered GRN
    filtered_edges = pd.DataFrame({
        'TF': [graph.vs[edge.source]['name'] for edge in graph.es],
        'target': [graph.vs[edge.target]['name'] for edge in graph.es],
        'importance': graph.es['weight'],
        'cluster': [graph.vs[edge.source]['cluster'] for edge in graph.es]
    })
    filtered_edges.to_csv(output_file, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Leiden clustering on the GRN")
    parser.add_argument("--grn_file", required=True, help="Path to the inferred GRN file")
    parser.add_argument("--resolution", type=float, required=True, help="Resolution parameter for Leiden clustering")
    parser.add_argument("--output_file", required=True, help="Path to save the filtered GRN")
    args = parser.parse_args()

    run_leiden_clustering(args.grn_file, args.resolution, args.output_file)
