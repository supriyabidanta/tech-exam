# preproc.py

import scanpy as sc
import logging

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[
                        logging.FileHandler("analysis.log", mode='w'),  # Overwrite the log file on each run
                        logging.StreamHandler()
                    ])

logger = logging.getLogger(__name__)


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
            sample_adata.var_names_make_unique()  # Ensure variable names are unique
            adatas[sample_id] = sample_adata

        adata = sc.concat(adatas, label="sample")
        adata.obs_names_make_unique()

        logger.info("Data loaded and variable names made unique for all samples.")

        return adata

class QualityControl:
    def filter_cells(self, adata, min_genes, max_genes):
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_cells(adata, max_genes=max_genes)
        logger.info("Cell filtering applied. Min genes: %s, Max genes: %s", min_genes, max_genes)

    def filter_genes(self, adata, min_cells):
        sc.pp.filter_genes(adata, min_cells=min_cells)
        logger.info("Gene filtering applied. Min cells per gene: %s", min_cells)

    def calculate_qc_metrics(self, adata, mt_prefix='mt'):
        mt_gene_mask = adata.var_names.str.startswith(mt_prefix)
        if any(mt_gene_mask):
            adata.obs['pct_counts_mt'] = (adata[:, mt_gene_mask].X.sum(axis=1) / adata.X.sum(axis=1)) * 100
            logger.info("Mitochondrial QC metrics calculated.")
        else:
            adata.obs['pct_counts_mt'] = 0.0
            logger.warning("No mitochondrial genes found. QC metrics for mitochondrial content set to 0.")

        qc_vars = [mt_prefix] if any(mt_gene_mask) else []
        sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, log1p=False, inplace=True)
        
    def run_qc_checks(self, adata, min_genes_per_cell, max_genes_per_cell,
                      min_cells_per_gene):
        self.filter_cells(adata, min_genes_per_cell, max_genes_per_cell)
        self.filter_genes(adata, min_cells_per_gene)
        self.calculate_qc_metrics(adata)
        logger.info("Quality control checks completed.")
        return adata


class Normalization:
    def normalize_total(self, adata, target_sum):
        sc.pp.normalize_total(adata, target_sum=target_sum)
        logger.info("Normalization by total counts with target sum: %s", target_sum)

    def log_transform(self, adata):
        sc.pp.log1p(adata)
        logger.info("Logarithmic transformation applied.")

    def normalize_data(self, adata, target_sum):
        self.normalize_total(adata, target_sum)
        self.log_transform(adata)
        logger.info("Data normalization completed.")
        return adata


class FeatureSelection:
    def highly_variable_genes(self, adata, min_mean, max_mean, min_disp):
        sc.pp.highly_variable_genes(adata, min_mean=min_mean,
                                    max_mean=max_mean, min_disp=min_disp)
        logger.info("Selection of highly variable genes completed.")

    def regress_out(self, adata, variables):
        sc.pp.regress_out(adata, variables)
        logger.info("Regression out of variables: %s", variables)

    def scale(self, adata, max_value):
        sc.pp.scale(adata, max_value=max_value)
        logger.info("Data scaling applied with max value: %s", max_value)

    def select_features(self, adata, min_mean, max_mean, min_disp,
                        regress_out_variables, scale_max_value):
        self.highly_variable_genes(adata, min_mean, max_mean, min_disp)
        self.regress_out(adata, regress_out_variables)
        self.scale(adata, scale_max_value)
        logger.info("Feature selection completed.")
        return adata


class DimensionalityReduction:
    def run_pca(self, adata, n_comps):
        sc.pp.pca(adata, n_comps=n_comps)
        logger.info("PCA done, check PCA data: %s", 'X_pca' in adata.obsm)
        return adata

    def run_umap(self, adata, min_dist, n_neighbors, n_pcs):
        if 'X_pca' in adata.obsm:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
            sc.tl.umap(adata, min_dist=min_dist)
            logger.info("Neighbors set for UMAP: %s", 'neighbors' in adata.uns)
        else:
            logger.warning("PCA results missing, cannot compute neighbors!")
        return adata

class Clustering:
    def run_leiden(self, adata, resolution, flavor='igraph', n_iterations=2, directed=False):
        if 'neighbors' in adata.uns:
            sc.tl.leiden(adata, resolution=resolution, flavor=flavor, n_iterations=n_iterations, directed=directed)
            logger.info("Leiden clustering executed successfully.")
        else:
            logger.error("Neighbors data missing, cannot perform Leiden clustering")

    def run_clustering(self, adata, leiden_resolution):

        logger.info("Attempting Leiden clustering, check for neighbors in .uns: %s", 'neighbors' in adata.uns)
        if 'neighbors' in adata.uns:
            self.run_leiden(adata, leiden_resolution)
            logger.info("Leiden clustering completed successfully.")
        else:
            logger.error("Neighbors data missing, cannot perform Leiden clustering")
        return adata
