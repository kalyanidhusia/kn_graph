import json
from pathlib import Path

out = Path("notebooks/01_hpxkg_walkthrough.ipynb")
out.parent.mkdir(parents=True, exist_ok=True)

def md(cell_id, text):
    return {
        "cell_type": "markdown",
        "id": cell_id,
        "metadata": {},
        "source": text
    }

def code(cell_id, text):
    return {
        "cell_type": "code",
        "execution_count": None,
        "id": cell_id,
        "metadata": {},
        "outputs": [],
        "source": text
    }

nb = {
    "metadata": {
        "kernelspec": {
            "display_name": "Python (hpxkg)",
            "language": "python",
            "name": "hpxkg"
        },
        "language_info": {
            "name": "python",
            "version": "3.11"
        },
        "title": "HPX Interactome Knowledge Graph Walkthrough",
        "authors": [
            {"name": "Kalyani Dhusia"}
        ]
    },
    "nbformat": 4,
    "nbformat_minor": 5,
    "cells": [
        md(
            "intro001",
            """# HPX Interactome Knowledge Graph

This notebook is a walkthrough of a reproducible proteomics-to-knowledge-graph workflow built from mzIdentML files associated with the PRIDE dataset **PXD060789**.

## Project goal
To identify proteins associated with HPX conditions, compare them against transferrin and control conditions, and organize the resulting proteins into an interpretable graph linking:

**Condition → Protein → Biological Process**

## Main result
The workflow identifies an HPX-enriched subnetwork centered on **LRP1**, with additional proteins linked to uptake, iron homeostasis, heme handling, complement, coagulation, protease regulation, and lipoprotein biology.
"""
        ),
        code(
            "load001",
            """from pathlib import Path
import pandas as pd

DATA = Path("../data/processed")
META = Path("../data/metadata")
FIGS = Path("../figures")

psm = pd.read_csv(DATA / "psm_table.csv")
protein_summary = pd.read_csv(DATA / "protein_evidence_summary.csv")
annotated = pd.read_csv(DATA / "protein_evidence_annotated.csv")
comparison = pd.read_csv(DATA / "protein_condition_comparison.csv")
top_hpx = pd.read_csv(DATA / "top_hpx_candidates.csv")
proc_ann = pd.read_csv(META / "protein_process_annotations.csv")
graph_table = pd.read_csv(DATA / "protein_process_graph_table.csv")
nodes = pd.read_csv(DATA / "process_graph_nodes.csv")
edges = pd.read_csv(DATA / "process_graph_edges.csv")

print("Loaded all major project tables.")
"""
        ),
        md(
            "mzid001",
            """## 1. Parsed mzID output

The project begins by parsing 12 mzIdentML files into peptide-spectrum match (PSM) and protein-level evidence tables.
"""
        ),
        code(
            "mzid002",
            """print("PSM rows:", len(psm))
print("Protein summary rows:", len(protein_summary))
print("Annotated protein rows:", len(annotated))
print("Unique accessions:", annotated["protein_accession_clean"].nunique())
print("Unique genes:", annotated["gene_symbol"].dropna().nunique())

psm.head()
"""
        ),
        code(
            "mzid003",
            """protein_summary.sort_values("n_peptides", ascending=False).head(15)
"""
        ),
        md(
            "candidates001",
            """## 2. HPX candidate proteins

Proteins were compared across three condition classes:
- **HPX**
- **Tf**
- **control**

The comparison table flags proteins that are present in HPX and absent from controls, and calculates a simple HPX-over-background score.
"""
        ),
        code(
            "candidates002",
            """top_hpx.head(20)
"""
        ),
        md(
            "candidates003",
            """### Interpretation

A few proteins stand out immediately:

- **LRP1** is the strongest HPX-specific anchor in this analysis.
- **TFRC** is shared between HPX and transferrin conditions, which is biologically sensible for uptake-related comparison.
- **CP, HP, HPR, SERPINA1, A2M** add extracellular and heme/iron-related context.
- **C1QB, C1QC, FCN3, F12** support complement and coagulation-linked modules.
"""
        ),
        md(
            "proc001",
            """## 3. Protein-to-process annotation

Top HPX-associated proteins were assigned small, interpretable biological process labels through a semi-manual curation workflow, producing a simple knowledge graph layer.
"""
        ),
        code(
            "proc002",
            """proc_ann.head(20)
"""
        ),
        code(
            "proc003",
            """proc_ann["process_group"].value_counts()
"""
        ),
        md(
            "proc004",
            """These process groups provide the backbone of the knowledge graph:

- uptake
- iron homeostasis
- heme handling
- complement
- coagulation
- protease regulation
- lipoprotein biology
"""
        ),
        md(
            "graph001",
            """## 4. Graph-ready tables

The final graph links condition classes to proteins and proteins to curated biological processes.
"""
        ),
        code(
            "graph002",
            """print("Node count:", len(nodes))
print("Edge count:", len(edges))
print("\\nNode types:")
print(nodes["node_type"].value_counts())

nodes.head(10)
"""
        ),
        code(
            "graph003",
            """edges.head(15)
"""
        ),
        md(
            "fig001",
            """## 5. Final figure

Below is the publication-style figure generated from the processed graph tables.
"""
        ),
        code(
            "fig002",
            """from IPython.display import Image, display
display(Image(filename=str(FIGS / "figure1_hpx_knowledge_graph.png")))
"""
        ),
        md(
            "focus001",
            """## 6. Key proteins summary
"""
        ),
        code(
            "focus002",
            """focus_genes = ["LRP1", "TFRC", "CP", "HP", "HPR", "A2M", "SERPINA1", "F12", "C1QB", "C1QC", "FCN3"]
comparison.loc[comparison["gene_symbol"].isin(focus_genes)].sort_values(
    ["candidate_HPX_over_background_score", "HPX"],
    ascending=[False, False]
)
"""
        ),
        md(
            "end001",
            """## Take-home message

This project demonstrates a reproducible workflow for:

- parsing mzIdentML proteomics outputs in Python
- summarizing protein-level evidence across conditions
- annotating proteins with UniProt-derived identifiers and names
- comparing HPX, transferrin, and control conditions
- building a small knowledge graph linking proteins to biological processes
- generating interactive and publication-style visual outputs

The main biological result is an **HPX-enriched subnetwork centered on LRP1**, with linked uptake, iron homeostasis, heme handling, complement, coagulation, protease regulation, and lipoprotein-associated proteins.
"""
        )
    ]
}

with open(out, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=2)

print(f"Wrote: {out}")
