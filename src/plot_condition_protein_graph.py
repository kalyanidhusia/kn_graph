#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import networkx as nx
import plotly.graph_objects as go

NODES = Path("data/processed/graph_nodes.csv")
EDGES = Path("data/processed/graph_edges.csv")
OUT_HTML = Path("results/condition_protein_graph.html")


def main():
    nodes = pd.read_csv(NODES)
    edges = pd.read_csv(EDGES)

    G = nx.Graph()

    for _, row in nodes.iterrows():
        G.add_node(
            row["id"],
            label=row["label"],
            node_type=row["node_type"]
        )

    for _, row in edges.iterrows():
        G.add_edge(
            row["source"],
            row["target"],
            weight=row["weight"],
            edge_type=row["edge_type"]
        )

    pos = nx.spring_layout(G, seed=42, k=0.9)

    edge_x = []
    edge_y = []
    for source, target in G.edges():
        x0, y0 = pos[source]
        x1, y1 = pos[target]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        line=dict(width=1, color="#888"),
        hoverinfo="none",
        mode="lines"
    )

    node_x = []
    node_y = []
    text = []
    color = []
    size = []

    for node, attrs in G.nodes(data=True):
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

        label = attrs.get("label", node)
        node_type = attrs.get("node_type", "Protein")
        text.append(f"{label}<br>{node_type}")

        if node_type == "Condition":
            color.append("#d62728")
            size.append(28)
        else:
            color.append("#1f77b4")
            size.append(14)

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=[attrs.get("label", n) for n, attrs in G.nodes(data=True)],
        textposition="top center",
        hovertext=text,
        hoverinfo="text",
        marker=dict(
            color=color,
            size=size,
            line=dict(width=1, color="black")
        )
    )

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        title="HPX condition–protein network",
        showlegend=False,
        template="plotly_white",
        margin=dict(l=20, r=20, t=50, b=20)
    )

    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(OUT_HTML)
    print(f"Wrote: {OUT_HTML}")


if __name__ == "__main__":
    main()
