import argparse

import numpy as np
import scanpy as sc
from anndata import AnnData


def get_data(path: str):
    adata = sc.read_10x_mtx(path, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    return adata

def prepare_data(adata: AnnData):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    for column in adata.obs.columns:
        try:
            adata.obs[f'{column}_log10'] = np.log10(adata.obs[column] + 1)
            adata.obs[f'{column}_log2'] = np.log2(adata.obs[column] + 1)
        except (ValueError, AttributeError, TypeError):
            pass
    return adata

def filter_data(adata: AnnData, min_genes: int, min_cells: int, mt_pct=5, min_counts=2500):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs.n_genes_by_counts < min_counts, :]
    adata = adata[adata.obs.pct_counts_mt < mt_pct, :]
    return adata

def prepare_for_dim_reduction(adata: AnnData):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    return adata


def dim_reduction(adata: AnnData):
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    return adata


def clustering(adata: AnnData, res: list):
    for resolution in res:
        sc.tl.leiden(adata, resolution = resolution, key_added = f"leiden_{resolution}")
    return adata


def get_markers(adata: AnnData, res: list):
    for resolution in res:
        sc.tl.rank_genes_groups(adata, groupby = f'leiden_{resolution}', method='wilcoxon',
                                key_added = f'{resolution}', pts = True)
    return adata


def save_adata(adata: AnnData, out_file: str):
    adata.write_h5ad(out_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Basic preprocessing: scanpy')
    parser.add_argument('--data', type=str, required=True,
                        help='Path to the CellRanger output')
    parser.add_argument('--out_file', type=str, required=True,
                        help='Path to the output h5ad file')
    parser.add_argument('--res', nargs='+', default=[1.0],
                        help='List of the resolutions for the clustering purposes')
    parser.add_argument('--min_genes', type=int, required=False, default=200)
    parser.add_argument('--min_cells', type=int, required=False, default=3)
    args = parser.parse_args()
    args.res = [float(x) for x in args.res]
    print(args)

    adata = get_data(args.data)
    adata = prepare_data(adata)
    adata = filter_data(adata, args.min_genes, args.min_cells)
    adata = prepare_for_dim_reduction(adata)
    adata = dim_reduction(adata)
    adata = clustering(adata, args.res)
    adata = get_markers(adata, args.res)
    save_adata(adata, args.out_file)
