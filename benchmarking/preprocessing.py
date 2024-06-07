import os
import scanpy as sc, anndata as ad
import harmonypy
import matplotlib as plt
import seaborn as sns



dsets = ["GSM5114461_S6_A11", "GSM5114464_S7_D20", "GSM5114474_M3_E7"]
adatas = {}
DATADIR = "/Users/swen/Github/CSE185_FinalProject"
dataset = sc.read_10x_mtx(DATADIR, prefix="GSM5114461_S6_A11_", cache=True)
for ds in dsets:
    print(ds)
    adatas[ds] = sc.read_10x_mtx(DATADIR, prefix=ds+"_", cache=True)
    adatas[ds].obs_names_make_unique()
combined = ad.concat(adatas, label="dsets")
combined.obs_names_make_unique()

sc.pp.filter_cells(combined, min_genes=200)
sc.pp.filter_cells(combined, min_counts=1000)
sc.pp.filter_genes(combined, min_cells=5)
sc.pp.filter_genes(combined, min_counts=15)

mt_genes = combined.var_names.str.startswith('MT-')
combined.obs['mito_counts'] = combined[:, mt_genes].X.sum(axis=1).A1
combined.obs['total_counts'] = combined.X.sum(axis=1).A1
combined.obs['percent_mito'] = (combined.obs['mito_counts'] / combined.obs['total_counts']) * 100

# Delete cells more than 25% of mitochondiral genes
adata_filt = combined[combined.obs['percent_mito'] < 25, :]

# Further filtering based on total counts and detected genes
adata_filt = adata_filt[(adata_filt.obs['total_counts'] > 1000) & (adata_filt.obs['n_genes'] > 200), :]

sc.pp.normalize_per_cell(adata_filt, counts_per_cell_after=1e4) # normalize to 10,000 reads/cell
sc.pp.log1p(adata_filt)

sc.pp.highly_variable_genes(adata_filt, batch_key='dsets', n_top_genes=500)

genes = ["GCG", "TTR",  "IAPP",  "GHRL", "PPY", "COL3A1",
    "CPA1", "CLPS", "REG1A", "CTRB1", "CTRB2", "PRSS2", "CPA2", "KRT19", "INS","SST","CELA3A", "VTCN1"]

adata_var = adata_filt[:, (adata_filt.var.index.isin(genes) | adata_filt.var["highly_variable"])]

sc.pp.pca(adata_var, n_comps=20)
sc.pl.pca(adata_var, color="dsets")

import scanpy.external as sce
sce.pp.harmony_integrate(adata_var, 'dsets', theta=2, nclust=50,  max_iter_harmony = 10,  max_iter_kmeans=10)
adata_var.obsm['X_pca'] = adata_var.obsm['X_pca_harmony']


adata_var.write("/Users/swen/Github/CSE185_FinalProject/harmony_integrated.h5ad")

# Plot PCA after Harmony integration
sc.pl.pca(adata_var, color="dsets")











