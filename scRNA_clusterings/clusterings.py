import scanpy as sc, anndata as ad
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
import sc3s
import leidenalg

# Clustering Methods
def kmeans_clustering(adata, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(adata.X)
    return kmeans.labels_

def hierarchical_clustering(adata, n_clusters):
    hierarchical = AgglomerativeClustering(n_clusters=n_clusters).fit(adata.X.toarray())
    return hierarchical.labels_

def leiden_clustering(adata):
    sc.tl.leiden(adata)
    return adata.obs['leiden'].values

def sc3s_clustering(adata, n_clusters):
    sc3s.tl.consensus(adata, n_clusters)
    return adata.obs[f'sc3s_{n_clusters}']

# Benchmarking with Silhouette Score
def benchmark_clustering(adata, labels):
    score = silhouette_score(adata.X, labels)
    return score

# Visualization
def visualize_clustering(adata, labels, title, method='umap'):
    adata.obs['cluster'] = labels
    if method == 'umap':
        sc.tl.umap(adata)
        sc.pl.umap(adata, color='cluster', title=title)
    elif method == 'tsne':
        sc.tl.tsne(adata)
        sc.pl.tsne(adata, color='cluster', title=title)