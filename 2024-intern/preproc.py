# preproc.py

import scanpy as sc


class DataLoader:
    def __init__(self, example_data):
        self.example_data = example_data

    def load_data(self):
        samples = {
            "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
        }
        adatas = {}

        for sample_id, filename in samples.items():
            path = self.example_data.fetch(filename)
            sample_adata = sc.read_10x_h5(path)
            sample_adata.var_names_make_unique()
            adatas[sample_id] = sample_adata

        adata = sc.concat(adatas, label="sample")
        adata.obs_names_make_unique()

        return adata


class QualityControl:
    def filter_cells(self, adata, min_genes, max_genes):
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_cells(adata, max_genes=max_genes)

    def filter_genes(self, adata, min_cells):
        sc.pp.filter_genes(adata, min_cells=min_cells)

    def calculate_qc_metrics(self, adata):
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                                   log1p=False, inplace=True)

    def run_qc_checks(self, adata, min_genes_per_cell, max_genes_per_cell,
                      min_cells_per_gene):
        self.filter_cells(adata, min_genes_per_cell, max_genes_per_cell)
        self.filter_genes(adata, min_cells_per_gene)
        self.calculate_qc_metrics(adata)
        return adata


class Normalization:
    def normalize_total(self, adata, target_sum):
        sc.pp.normalize_total(adata, target_sum=target_sum)

    def log_transform(self, adata):
        sc.pp.log1p(adata)

    def normalize_data(self, adata, target_sum):
        self.normalize_total(adata, target_sum)
        self.log_transform(adata)
        return adata


class FeatureSelection:
    def highly_variable_genes(self, adata, min_mean, max_mean, min_disp):
        sc.pp.highly_variable_genes(adata, min_mean=min_mean,
                                    max_mean=max_mean, min_disp=min_disp)

    def regress_out(self, adata, variables):
        sc.pp.regress_out(adata, variables)

    def scale(self, adata, max_value):
        sc.pp.scale(adata, max_value=max_value)

    def select_features(self, adata, min_mean, max_mean, min_disp,
                        regress_out_variables, scale_max_value):
        self.highly_variable_genes(adata, min_mean, max_mean, min_disp)
        self.regress_out(adata, regress_out_variables)
        self.scale(adata, scale_max_value)
        return adata


class DimensionalityReduction:
    def run_pca(self, adata, n_comps):
        sc.pp.pca(adata, svd_solver='arpack', n_comps=n_comps)
        return adata

    def run_umap(self, adata, min_dist, n_neighbors, n_pcs):
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(adata, min_dist=min_dist)
        return adata


class Clustering:
    def run_leiden(self, adata, resolution):
        sc.tl.leiden(adata, resolution=resolution)

    def run_clustering(self, adata, leiden_resolution):
        self.run_leiden(adata, leiden_resolution)
        return adata
