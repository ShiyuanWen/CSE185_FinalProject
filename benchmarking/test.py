from scRNA_clusterings import (
    kmeans_clustering, hierarchical_clustering, leiden_clustering,
    sc3s_clustering, benchmark_clustering, visualize_clustering
)
import scanpy as sc
import time
start_time = time.time()

adata_var = sc.read("harmony_integrated.h5ad")
sc.pp.neighbors(adata_var)


kmeans_labels = leiden_clustering(adata_var)
kmeans_score = benchmark_clustering(adata_var, kmeans_labels)


visualize_clustering(adata_var, kmeans_labels, 'Leiden', method='umap')

print(f'sc3s-clustering score Silhouette Score: {kmeans_score}')

print(adata_var.obs.columns)

sc.tl.rank_genes_groups(adata_var, groupby='leiden', method='t-test')
rank_genes_groups_df = sc.get.rank_genes_groups_df(adata_var, group=None)

rank_genes_groups_df = sc.get.rank_genes_groups_df(adata_var, group='1')

top_10_genes_group_5 = rank_genes_groups_df.head(10)

print("Top 10 genes for group 6:")
print(top_10_genes_group_5[['names', 'logfoldchanges', 'pvals_adj']])

end_time = time.time()

execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")
