import scanpy as sc, anndata as ad
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
import sc3s
import leidenalg
import scipy.sparse


# Clustering Methods
def kmeans_clustering(adata, n_clusters):
    # Perform kmeans clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(adata.X)
    # Add clustering labels into adata object
    adata.obs['kmeans'] = kmeans.labels_.astype(str)
    adata.obs['kmeans'] = adata.obs['kmeans'].astype('category')
    return adata.obs['kmeans'].values

def hierarchical_clustering(adata, n_clusters):
    # Perform hierarchical clustering
    if scipy.sparse.issparse(adata.X):
        data = adata.X.toarray()
    else:
        data = adata.X
    hierarchical = AgglomerativeClustering(n_clusters=n_clusters).fit(data)
    # Add clustering labels into adata object
    adata.obs['hierarchical'] = hierarchical.labels_.astype(str)
    adata.obs['hierarchical'] = adata.obs['hierarchical'].astype('category')
    return adata.obs['hierarchical'].values

def leiden_clustering(adata):
    # Perform leiden clustering
    sc.tl.leiden(adata)
    # Add clustering label into adata object
    return adata.obs['leiden'].values

def sc3s_clustering(adata, n_clusters):
    # Perform sc3s_clustering
    sc3s.tl.consensus(adata, n_clusters)
    # Add clustering label into adata object
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