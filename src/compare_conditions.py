#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

INFILE = Path("data/processed/protein_evidence_with_metadata.csv")
OUTDIR = Path("data/processed")
OUTDIR.mkdir(parents=True, exist_ok=True)


def main():
    if not INFILE.exists():
        raise FileNotFoundError(f"Missing input file: {INFILE}")

    df = pd.read_csv(INFILE)

    # Keep only non-contaminants
    df = df.loc[~df["is_contaminant"].fillna(False)].copy()

    # Require at least 2 unique peptides for a more stable first-pass set
    df = df.loc[df["n_peptides"] >= 2].copy()

    summary = (
        df.groupby(
            ["protein_accession_clean", "gene_symbol", "protein_name_final", "condition_class"],
            dropna=False,
            as_index=False
        )
        .agg(
            files_detected=("source_file", "nunique"),
            total_psms=("n_psms", "sum"),
            max_peptides=("n_peptides", "max")
        )
    )

    pivot = (
        summary.pivot_table(
            index=["protein_accession_clean", "gene_symbol", "protein_name_final"],
            columns="condition_class",
            values="files_detected",
            fill_value=0
        )
        .reset_index()
    )

    for col in ["HPX", "Tf", "control"]:
        if col not in pivot.columns:
            pivot[col] = 0

    pivot["present_in_HPX"] = pivot["HPX"] > 0
    pivot["present_in_Tf"] = pivot["Tf"] > 0
    pivot["present_in_control"] = pivot["control"] > 0

    pivot["candidate_HPX_specific"] = pivot["present_in_HPX"] & (~pivot["present_in_control"])
    pivot["candidate_shared_HPX_Tf"] = pivot["present_in_HPX"] & pivot["present_in_Tf"]
    pivot["candidate_HPX_over_background_score"] = pivot["HPX"] - pivot["control"]

    pivot = pivot.sort_values(
        ["candidate_HPX_specific", "candidate_HPX_over_background_score", "HPX"],
        ascending=[False, False, False]
    )

    comparison_file = OUTDIR / "protein_condition_comparison.csv"
    top_file = OUTDIR / "top_hpx_candidates.csv"

    pivot.to_csv(comparison_file, index=False)
    pivot.loc[pivot["present_in_HPX"]].head(100).to_csv(top_file, index=False)

    print(f"Wrote: {comparison_file}")
    print(f"Wrote: {top_file}")

    print("\nTop HPX candidates:")
    show_cols = [
        "protein_accession_clean",
        "gene_symbol",
        "protein_name_final",
        "HPX",
        "Tf",
        "control",
        "candidate_HPX_specific",
        "candidate_shared_HPX_Tf",
        "candidate_HPX_over_background_score"
    ]
    print(
        pivot.loc[pivot["present_in_HPX"], show_cols]
        .head(20)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
