#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

ANNOTATED = Path("data/processed/protein_evidence_annotated.csv")
METADATA = Path("data/metadata/sample_metadata.csv")
OUTFILE = Path("data/processed/protein_evidence_with_metadata.csv")


def main():
    if not ANNOTATED.exists():
        raise FileNotFoundError(f"Missing annotated file: {ANNOTATED}")
    if not METADATA.exists():
        raise FileNotFoundError(f"Missing metadata file: {METADATA}")

    prot = pd.read_csv(ANNOTATED)
    meta = pd.read_csv(METADATA)

    merged = prot.merge(meta, on="source_file", how="left")

    print("Protein rows:", len(prot))
    print("Metadata rows:", len(meta))
    print("Merged rows:", len(merged))
    print("Rows missing metadata:", merged["condition_short"].isna().sum())

    merged.to_csv(OUTFILE, index=False)
    print(f"Wrote: {OUTFILE}")

    print("\nCondition counts:")
    print(merged["condition_short"].value_counts(dropna=False))

    print("\nCondition class counts:")
    print(merged["condition_class"].value_counts(dropna=False))

    print("\nTop rows:")
    show_cols = [
        "source_file",
        "protein_accession",
        "gene_symbol",
        "protein_name_final",
        "n_peptides",
        "condition_short",
        "condition_class",
        "batch",
        "is_contaminant",
    ]
    show_cols = [c for c in show_cols if c in merged.columns]
    print(merged[show_cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
