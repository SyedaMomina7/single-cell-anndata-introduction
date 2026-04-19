# Single-Cell Data Management using AnnData
### Getting Started Tutorial

---

## Overview

This repository demonstrates a complete workflow for working with **AnnData**,
a core data structure used in **single-cell RNA sequencing (scRNA-seq)** analysis in Python.

AnnData is widely used in computational biology for managing high-dimensional
single-cell datasets efficiently. It is the backbone of the **Scanpy ecosystem**.

---

## Objectives

- Understand how AnnData structures single-cell datasets
- Store and manage gene expression matrices
- Attach metadata to cells and genes
- Subset and filter data efficiently
- Store embeddings (PCA, UMAP, etc.)
- Manage multiple data layers
- Save and reload `.h5ad` files

---

## AnnData Structure Overview
AnnData object with n_obs × n_vars = 100 × 2000

| Component | Shape             | Description                          |
|-----------|-------------------|--------------------------------------|
| `.X`      | `n_obs × n_vars`  | Primary gene expression matrix       |
| `.obs`    | `n_obs × ...`     | Cell-level metadata (DataFrame)      |
| `.var`    | `n_vars × ...`    | Gene-level metadata (DataFrame)      |
| `.obsm`   | `n_obs × k`       | Embeddings (e.g., PCA, UMAP)         |
| `.varm`   | `n_vars × k`      | Gene-level annotations or embeddings |
| `.uns`    | unstructured      | Unstructured metadata (dict)         |
| `.layers` | `n_obs × n_vars`  | Alternative data representations     |

---

## Step-by-Step Workflow

### 1. Import Libraries

```python
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
```

---

### 2. Create an AnnData Object

Simulate a sparse gene expression matrix using a Poisson distribution:

```python
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)))
adata = ad.AnnData(counts)
```

- **100** observations (cells) — rows
- **2000** variables (genes) — columns
- Stored as a **sparse matrix** (`csr_matrix`) for memory efficiency

---

### 3. Inspect the AnnData Object

```python
adata
# AnnData object with n_obs × n_vars = 100 × 2000
```

This displays the shape, layers, and metadata slots currently populated.

---

### 4. Add Cell Metadata (`obs`)

```python
adata.obs["cell_type"] = np.random.choice(
    ["B cell", "T cell", "Monocyte"],
    size=adata.n_obs
)
```

> `.obs` is a **pandas DataFrame** indexed by cell barcodes.  
> Each row corresponds to one cell.

---

### 5. Add Gene Metadata (`var`)

```python
adata.var["gene_id"] = [f"Gene_{i}" for i in range(adata.n_vars)]
```

> `.var` is also a **pandas DataFrame** indexed by gene names.  
> Each row corresponds to one gene.

---

### 6. Subset Cells

```python
adata_bcells = adata[adata.obs.cell_type == "B cell"].copy()
```

> Subsetting returns a **view** of the original object.  
> `.copy()` creates an independent object to avoid modifying the original.

---

### 7. Store Embeddings (`obsm`)

```python
adata.obsm["X_pca"] = np.random.normal(size=(adata.n_obs, 50))
adata.obsm["X_umap"] = np.random.normal(size=(adata.n_obs, 2))
```

> `.obsm` stores **multi-dimensional arrays** aligned to cells.  
> Convention: keys are prefixed with `X_` (e.g., `X_umap`, `X_pca`).

---

### 8. Store Gene-Level Features (`varm`)

```python
adata.varm["gene_embedding"] = np.random.normal(size=(adata.n_vars, 5))
```

> `.varm` stores multi-dimensional arrays aligned to genes.

---

### 9. Store Unstructured Metadata (`uns`)

```python
adata.uns["dataset_info"] = {
    "source": "AnnData tutorial",
    "type": "synthetic dataset",
    "anndata_version": ad.__version__
}
```

> `.uns` is a free-form **dictionary** for any metadata that doesn't
> fit into `.obs` or `.var` — e.g., color palettes, pipeline parameters.

---

### 10. Store Alternative Data Representations (`layers`)

```python
import scipy.sparse as sp

adata.layers["counts"] = adata.X.copy()           # raw counts
adata.layers["log1p"] = np.log1p(adata.X.toarray())  # log-normalized
```

> `.layers` allows storing **multiple versions** of the expression matrix
> (raw counts, normalized, log-transformed) all in one object.

---

### 11. Convert to DataFrame (Optional)

```python
df = adata.to_df()
df.head()
```

> Converts `.X` to a dense **pandas DataFrame** with cell and gene labels.
> Use only for small subsets — large matrices may exceed memory.

---

### 12. Save AnnData Object

```python
adata.write_h5ad("results/results.h5ad")
```

> Saves the entire AnnData object to disk in **HDF5 format** (`.h5ad`),
> including the expression matrix, all metadata, embeddings, and layers.

> ⚠️ Note: `ad.read()` is deprecated. Use `ad.write_h5ad()` to save.

---

### 13. Load Saved Object

```python
adata = ad.read_h5ad("results/results.h5ad")
```

> ⚠️ Note: `ad.read()` is deprecated as of AnnData 0.10+.  
> Use `ad.read_h5ad()` instead.

---
## Results  

- Successfully created an AnnData object with a sparse matrix  
- Assigned meaningful observation and variable names  
- Added metadata using `obs` and `var`  

- Performed subsetting using:
  - Index-based selection  
  - Metadata filtering  

- Stored multi-dimensional data in:
  - `obsm` (e.g., UMAP embeddings)  
  - `varm`  

- Used `uns` for unstructured metadata storage  

- Implemented layers to store transformed data (e.g., log-transformed values)  

- Demonstrated the concept of **views vs copies** and memory-efficient operations  

- Saved and loaded data using `.h5ad` format  
## Key Concepts

### Why AnnData is Powerful

| Feature | Benefit |
|--------|---------|
| Sparse matrix support | Handles millions of cells efficiently |
| Aligned metadata | `.obs` and `.var` stay in sync with `.X` |
| Multiple layers | Store raw, normalized, and transformed data together |
| Scanpy integration | Drop-in compatibility with the full sc analysis ecosystem |
| `.h5ad` format | Portable, compressed, and fast I/O |


---
## Workflow Summary

Raw Count Matrix  
↓  
AnnData Object (.X)  
↓  
Add Metadata (.obs / .var)  
↓  
Subset / Filter Cells  
↓  
Embeddings (.obsm: PCA / UMAP)  
↓  
Layers (.layers: log / normalized)  
↓  
Save to Disk (.h5ad)

## Repository Structure

.
├── ann data project/
│   └── anndata_getting_started.ipynb
├── README.md
└── requirements.txt


## Requirements

```txt
anndata>=0.10.0
numpy
pandas
scipy
```

Install with:

```bash
pip install anndata numpy pandas scipy
```

---

## Final Notes

This tutorial covers the **foundational structure** of single-cell data analysis
using AnnData. It is a core building block for more advanced workflows:

- Quality control (QC) filtering
- Normalization and scaling
- PCA / UMAP visualization
- Clustering (Leiden / Louvain algorithms)
- Differential expression analysis

---

## Future Improvements

- [ ] Integrate full Scanpy pipeline (`sc.pp`, `sc.tl`, `sc.pl`)
- [ ] Add real scRNA-seq dataset (e.g., PBMC 3k from 10x Genomics)
- [ ] Include UMAP and clustering plots
- [ ] Add reproducible environment (`conda` or `Docker`)
- [ ] Migrate to `AnnData` 0.13+ features (e.g., `ExperimentalAnnCollection`)

---

##Author

Syeda Momina Assad

## References

- [AnnData Documentation](https://anndata.readthedocs.io)

