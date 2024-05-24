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
kmeans_labels = kmeans_clustering(adata_var, n_clusters=7)
kmeans_score = benchmark_clustering(adata_var, kmeans_labels)
print(f'K-Means Silhouette Score: {kmeans_score}')

visualize_clustering(adata_var, kmeans_labels, 'K-Means Clustering', method='umap')
