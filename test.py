from clusterings import (
    kmeans_clustering, hierarchical_clustering, leiden_clustering,
    sc3s_clustering, benchmark_clustering, visualize_clustering
)
import scanpy as sc, anndata as ad
import matplotlib as plt
import seaborn as sns
import leidenalg

adata_var = sc.read("/Users/swen/Github/CSE185_FinalProject/harmony_integrated.h5ad")
sc.pp.neighbors(adata_var)


sc3s_labels = sc3s_clustering(adata_var, n_clusters=8)
sc3s_score = benchmark_clustering(adata_var, sc3s_labels)


visualize_clustering(adata_var, sc3s_labels, 'sc3s', method='umap')

print(f'sc3s-clustering score Silhouette Score: {sc3s_score}')
