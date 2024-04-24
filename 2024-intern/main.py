# main.py

import argparse
import pooch
from preproc import (
    DataLoader,
    QualityControl,
    Normalization,
    FeatureSelection,
    DimensionalityReduction,
    Clustering
)


def main(output_file):
    example_data = pooch.create(
        path=pooch.os_cache("scverse_tutorials"),
        base_url="doi:10.6084/m9.figshare.22716739.v1/",
    )
    example_data.load_registry_from_doi()

    min_genes_per_cell = 200
    max_genes_per_cell = 2500
    min_cells_per_gene = 3

    target_sum = 1e4

    min_mean = 0.0125
    max_mean = 3
    min_disp = 0.5
    regress_out_variables = ['total_counts', 'pct_counts_mt']
    scale_max_value = 10

    n_pcs = 300
    n_neighbors = 15
    umap_min_dist = 0.5

    leiden_resolution = 2.3

    dl = DataLoader(example_data)
    adata = dl.load_data()

    qc = QualityControl()
    adata = qc.run_qc_checks(adata, min_genes_per_cell, max_genes_per_cell,
                             min_cells_per_gene)

    normalization = Normalization()
    adata = normalization.normalize_data(adata, target_sum)

    feature_selection = FeatureSelection()
    adata = feature_selection.select_features(
        adata,
        min_mean,
        max_mean,
        min_disp,
        regress_out_variables,
        scale_max_value
    )

    clustering = Clustering()
    adata = clustering.run_clustering(adata, leiden_resolution)

    dr = DimensionalityReduction()
    adata = dr.run_pca(adata, n_pcs)
    adata = dr.run_umap(adata, umap_min_dist, n_neighbors, n_pcs)

    adata.write(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Preprocess single-cell RNA-seq data.'
    )
    parser.add_argument('--output_file', type=str, required=True,
                        help='Path to save the preprocessed data.')

    args = parser.parse_args()

    main(args.output_file)

