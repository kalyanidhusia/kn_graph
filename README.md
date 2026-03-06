# HPX Interactome Knowledge Graph

A reproducible proteomics-to-graph workflow built from PRIDE-linked mzIdentML files and supplementary study tables for the HPX interactome project (PXD060789).

## Goal
Parse mzIdentML files, summarize peptide/protein evidence, integrate study metadata and supplementary biological annotations, and generate graph-ready tables for knowledge graph visualization of HPX-associated interactomes.

## Current workflow
1. Parse `.mzid` files into peptide-spectrum/protein evidence tables
2. Summarize protein-level evidence by condition
3. Merge with supplementary biological annotations
4. Build node/edge tables for graph analysis
5. Visualize interactome modules related to iron homeostasis, receptor-mediated endocytosis, immunity, and coagulation

## Repository structure
- `data/raw/mzid/` raw mzIdentML files
- `data/metadata/` sample metadata and annotation tables
- `data/processed/` parsed outputs
- `notebooks/` exploratory notebooks
- `src/` reusable parsing and processing scripts
- `results/` outputs for figures and graph files

## Environment
```bash
conda env create -f environment.yml
conda activate hpxkg
python -m ipykernel install --user --name hpxkg --display-name "Python (hpxkg)"
