#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

INFILE = Path("data/processed/protein_condition_comparison.csv")
OUTDIR = Path("data/processed")
OUTDIR.mkdir(parents=True, exist_ok=True)

TOP_N = 50


def main():
    df = pd.read_csv(INFILE)

    # Focus on graph-worthy proteins
    graph_df = df.loc[
        (df["present_in_HPX"] == True) &
        (df["present_in_control"] == False)
    ].copy()

    graph_df = graph_df.sort_values(
        ["candidate_HPX_over_background_score", "HPX", "Tf"],
        ascending=[False, False, False]
    ).head(TOP_N)

    # Nodes
    protein_nodes = graph_df[[
        "protein_accession_clean", "gene_symbol", "protein_name_final"
    ]].copy()
    protein_nodes = protein_nodes.rename(columns={
        "protein_accession_clean": "id"
    })
    protein_nodes["label"] = protein_nodes["gene_symbol"].fillna(protein_nodes["id"])
    protein_nodes["node_type"] = "Protein"

    condition_nodes = pd.DataFrame({
        "id": ["HPX", "Tf", "control"],
        "label": ["HPX", "Tf", "control"],
        "gene_symbol": [None, None, None],
        "protein_name_final": [None, None, None],
        "node_type": ["Condition", "Condition", "Condition"]
    })

    nodes = pd.concat([
        condition_nodes[["id", "label", "gene_symbol", "protein_name_final", "node_type"]],
        protein_nodes[["id", "label", "gene_symbol", "protein_name_final", "node_type"]],
    ], ignore_index=True)

    # Edges
    edge_rows = []
    for _, row in graph_df.iterrows():
        protein_id = row["protein_accession_clean"]

        if row.get("HPX", 0) > 0:
            edge_rows.append({
                "source": "HPX",
                "target": protein_id,
                "edge_type": "detected_in",
                "weight": row["HPX"]
            })

        if row.get("Tf", 0) > 0:
            edge_rows.append({
                "source": "Tf",
                "target": protein_id,
                "edge_type": "detected_in",
                "weight": row["Tf"]
            })

        if row.get("control", 0) > 0:
            edge_rows.append({
                "source": "control",
                "target": protein_id,
                "edge_type": "detected_in",
                "weight": row["control"]
            })

    edges = pd.DataFrame(edge_rows)

    nodes.to_csv(OUTDIR / "graph_nodes.csv", index=False)
    edges.to_csv(OUTDIR / "graph_edges.csv", index=False)

    print("Wrote:", OUTDIR / "graph_nodes.csv")
    print("Wrote:", OUTDIR / "graph_edges.csv")
    print("\nTop graph proteins:")
    print(graph_df[[
        "protein_accession_clean",
        "gene_symbol",
        "protein_name_final",
        "HPX",
        "Tf",
        "control",
        "candidate_HPX_specific",
        "candidate_shared_HPX_Tf"
    ]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
