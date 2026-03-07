#!/usr/bin/env python3

from pathlib import Path
import textwrap

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


GRAPH_TABLE = Path("data/processed/protein_process_graph_table.csv")
COMPARISON = Path("data/processed/protein_condition_comparison.csv")
OUTDIR = Path("figures")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_PNG = OUTDIR / "figure1_hpx_knowledge_graph.png"
OUT_PDF = OUTDIR / "figure1_hpx_knowledge_graph.pdf"

TOP_N_BAR = 12


def wrap_label(s, width=22):
    if pd.isna(s):
        return ""
    return "\n".join(textwrap.wrap(str(s), width=width))


def prepare_graph_df():
    df = pd.read_csv(GRAPH_TABLE)

    # Prune background/noisy groups for a cleaner publication figure
    drop_groups = {"immune_background", "metabolism"}
    df = df.loc[~df["process_group"].isin(drop_groups)].copy()

    # Keep a focused protein set:
    # HPX present, control absent
    df = df.loc[
        (df["HPX"] > 0) &
        (df["control"] == 0)
    ].copy()

    # Optional: sort for deterministic behavior
    df = df.sort_values(
        ["HPX", "Tf", "gene_symbol", "process_group", "process_label"],
        ascending=[False, False, True, True, True]
    )

    return df


def build_graph(graph_df):
    G = nx.Graph()

    # condition nodes
    for cond in ["HPX", "Tf"]:
        G.add_node(cond, label=cond, node_type="Condition")

    # protein nodes
    for _, row in graph_df[[
        "protein_accession_clean", "gene_symbol"
    ]].drop_duplicates().iterrows():
        pid = row["protein_accession_clean"]
        label = row["gene_symbol"] if pd.notna(row["gene_symbol"]) else pid
        G.add_node(pid, label=label, node_type="Protein")

    # process nodes
    for _, row in graph_df[[
        "process_label", "process_group"
    ]].drop_duplicates().iterrows():
        proc = row["process_label"]
        G.add_node(proc, label=proc, node_type="Process", process_group=row["process_group"])

    # condition -> protein edges
    for _, row in graph_df[[
        "protein_accession_clean", "HPX", "Tf"
    ]].drop_duplicates().iterrows():
        pid = row["protein_accession_clean"]

        if row["HPX"] > 0:
            G.add_edge("HPX", pid, edge_type="detected_in", weight=float(row["HPX"]))

        if row["Tf"] > 0:
            G.add_edge("Tf", pid, edge_type="detected_in", weight=float(row["Tf"]))

    # protein -> process edges
    for _, row in graph_df[[
        "protein_accession_clean", "process_label"
    ]].drop_duplicates().iterrows():
        G.add_edge(
            row["protein_accession_clean"],
            row["process_label"],
            edge_type="involved_in",
            weight=1.0
        )

    return G


def draw_graph_panel(ax, graph_df):
    G = build_graph(graph_df)

    pos = nx.spring_layout(G, seed=42, k=1.15)

    # Node groups
    condition_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "Condition"]
    protein_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "Protein"]
    process_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "Process"]

    # Edge groups
    cond_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("edge_type") == "detected_in"]
    proc_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("edge_type") == "involved_in"]

    # Draw edges
    nx.draw_networkx_edges(
        G, pos, edgelist=cond_edges, ax=ax,
        width=1.8, alpha=0.65, edge_color="#9aa0a6"
    )
    nx.draw_networkx_edges(
        G, pos, edgelist=proc_edges, ax=ax,
        width=1.2, alpha=0.5, edge_color="#c7c7c7", style="dashed"
    )

    # Process node colors by group
    process_color_map = {
        "uptake": "#2ca02c",
        "iron_homeostasis": "#ff7f0e",
        "heme_handling": "#d62728",
        "complement": "#9467bd",
        "coagulation": "#8c564b",
        "protease_regulation": "#17becf",
        "lipoprotein_biology": "#e377c2",
    }

    proc_colors = []
    for n in process_nodes:
        group = G.nodes[n].get("process_group")
        proc_colors.append(process_color_map.get(group, "#7f7f7f"))

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos, nodelist=condition_nodes, ax=ax,
        node_color="#111111", node_size=1800, linewidths=1.0, edgecolors="black"
    )
    nx.draw_networkx_nodes(
        G, pos, nodelist=protein_nodes, ax=ax,
        node_color="#4c78a8", node_size=900, linewidths=0.8, edgecolors="black"
    )
    nx.draw_networkx_nodes(
        G, pos, nodelist=process_nodes, ax=ax,
        node_color=proc_colors, node_size=1300, linewidths=0.8, edgecolors="black"
    )

    # Labels
    labels = {}
    for n, d in G.nodes(data=True):
        label = d.get("label", n)
        if d.get("node_type") == "Process":
            label = wrap_label(label, width=20)
        labels[n] = label

    nx.draw_networkx_labels(
        G, pos, labels=labels, ax=ax,
        font_size=8, font_weight="bold"
    )

    ax.set_title("A. Condition → Protein → Process knowledge graph", loc="left", fontweight="bold")
    ax.axis("off")


def draw_bar_panel(ax):
    df = pd.read_csv(COMPARISON)

    # Focused and interpretable subset
    df = df.loc[
        (df["present_in_HPX"] == True) &
        (df["present_in_control"] == False)
    ].copy()

    # Remove likely noisy background from the first figure
    drop_genes = {"IGKC", "IGHG3", "IGHA1", "IGHA2", "IGLC6"}
    df = df.loc[~df["gene_symbol"].isin(drop_genes)].copy()

    df = df.sort_values(
        ["candidate_HPX_over_background_score", "HPX", "Tf"],
        ascending=[False, False, False]
    ).head(TOP_N_BAR)

    labels = df["gene_symbol"].fillna(df["protein_accession_clean"]).tolist()
    scores = df["candidate_HPX_over_background_score"].tolist()
    colors = ["#d62728" if x == "LRP1" else "#4c78a8" for x in labels]

    ax.barh(labels[::-1], scores[::-1], color=colors[::-1], edgecolor="black", linewidth=0.6)
    ax.set_xlabel("HPX over background score")
    ax.set_title("B. Top HPX-enriched proteins", loc="left", fontweight="bold")

    # Annotate bars with HPX/Tf counts
    for i, (_, row) in enumerate(df.iloc[::-1].iterrows()):
        txt = f'HPX={int(row["HPX"])}, Tf={int(row["Tf"])}'
        ax.text(
            row["candidate_HPX_over_background_score"] + 0.03,
            i,
            txt,
            va="center",
            fontsize=8
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    graph_df = prepare_graph_df()

    fig, axes = plt.subplots(
        1, 2, figsize=(16, 8),
        gridspec_kw={"width_ratios": [1.5, 1.0]}
    )

    draw_graph_panel(axes[0], graph_df)
    draw_bar_panel(axes[1])

    fig.suptitle(
        "HPX-associated proteomics knowledge graph identifies LRP1-centered uptake and extracellular functional modules",
        fontsize=14,
        fontweight="bold",
        y=0.98
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"Wrote: {OUT_PNG}")

    try:
        fig.savefig(OUT_PDF, bbox_inches="tight")
        print(f"Wrote: {OUT_PDF}")
    except Exception as e:
        print(f"PDF save skipped: {e}")

    plt.close(fig)


if __name__ == "__main__":
    main()
