#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

COMPARISON = Path("data/processed/protein_condition_comparison.csv")
ANNOTATIONS = Path("data/metadata/protein_process_annotations.csv")
OUTDIR = Path("data/processed")
OUTDIR.mkdir(parents=True, exist_ok=True)


def main():
    if not COMPARISON.exists():
        raise FileNotFoundError(f"Missing file: {COMPARISON}")
    if not ANNOTATIONS.exists():
        raise FileNotFoundError(f"Missing file: {ANNOTATIONS}")

    comp = pd.read_csv(COMPARISON)
    ann = pd.read_csv(ANNOTATIONS)

    # Focus on proteins present in HPX and absent from controls
    graph_proteins = comp.loc[
        (comp["present_in_HPX"] == True) &
        (comp["present_in_control"] == False)
    ].copy()

    graph_df = graph_proteins.merge(
        ann,
        on=["protein_accession_clean", "gene_symbol", "protein_name_final"],
        how="left"
    )

    # Keep only annotated proteins for the process graph
    graph_df = graph_df.loc[graph_df["process_label"].notna()].copy()

    # ---------- Nodes ----------
    condition_nodes = pd.DataFrame({
        "id": ["HPX", "Tf", "control"],
        "label": ["HPX", "Tf", "control"],
        "node_type": ["Condition", "Condition", "Condition"]
    })

    protein_nodes = (
        graph_df[["protein_accession_clean", "gene_symbol", "protein_name_final"]]
        .drop_duplicates()
        .rename(columns={"protein_accession_clean": "id"})
    )
    protein_nodes["label"] = protein_nodes["gene_symbol"].fillna(protein_nodes["id"])
    protein_nodes["node_type"] = "Protein"

    process_nodes = (
        graph_df[["process_label", "process_group"]]
        .drop_duplicates()
        .rename(columns={"process_label": "id"})
    )
    process_nodes["label"] = process_nodes["id"]
    process_nodes["node_type"] = "Process"

    nodes = pd.concat([
        condition_nodes[["id", "label", "node_type"]],
        protein_nodes[["id", "label", "node_type"]],
        process_nodes[["id", "label", "node_type"]],
    ], ignore_index=True).drop_duplicates()

    # ---------- Edges ----------
    edge_rows = []

    # Condition -> Protein
    for _, row in graph_df.drop_duplicates(
        subset=[
            "protein_accession_clean", "gene_symbol", "protein_name_final",
            "HPX", "Tf", "control"
        ]
    ).iterrows():
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

    # Protein -> Process
    for _, row in graph_df.drop_duplicates(
        subset=["protein_accession_clean", "process_label"]
    ).iterrows():
        edge_rows.append({
            "source": row["protein_accession_clean"],
            "target": row["process_label"],
            "edge_type": "involved_in",
            "weight": 1
        })

    edges = pd.DataFrame(edge_rows).drop_duplicates()

    nodes_file = OUTDIR / "process_graph_nodes.csv"
    edges_file = OUTDIR / "process_graph_edges.csv"
    graph_file = OUTDIR / "protein_process_graph_table.csv"

    nodes.to_csv(nodes_file, index=False)
    edges.to_csv(edges_file, index=False)
    graph_df.to_csv(graph_file, index=False)

    print(f"Wrote: {nodes_file}")
    print(f"Wrote: {edges_file}")
    print(f"Wrote: {graph_file}")

    print("\nAnnotated graph proteins:")
    print(
        graph_df[
            [
                "protein_accession_clean",
                "gene_symbol",
                "protein_name_final",
                "HPX",
                "Tf",
                "control",
                "process_label",
                "process_group"
            ]
        ]
        .drop_duplicates()
        .sort_values(["HPX", "Tf", "gene_symbol"], ascending=[False, False, True])
        .head(30)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
