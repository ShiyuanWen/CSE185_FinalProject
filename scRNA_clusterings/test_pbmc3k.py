import scanpy as sc
from scRNA_clusterings import (
    kmeans_clustering, hierarchical_clustering, leiden_clustering,
    sc3s_clustering, benchmark_clustering, visualize_clustering
)

adata = sc.datasets.pbmc3k()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden', title="leiden")

sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
rank_genes_groups_df = sc.get.rank_genes_groups_df(adata, group=None)

rank_genes_groups_df = sc.get.rank_genes_groups_df(adata, group='4')

top_10_genes_group_5 = rank_genes_groups_df.head(10)


print("Top 10 genes for group 6:")
print(top_10_genes_group_5[['names', 'pvals_adj']])


