#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

INFILE = Path("data/processed/top_hpx_candidates.csv")
OUTFILE = Path("data/metadata/protein_process_annotations_template.csv")

TOP_N = 30


def main():
    if not INFILE.exists():
        raise FileNotFoundError(f"Missing input file: {INFILE}")

    df = pd.read_csv(INFILE).head(TOP_N).copy()

    template = df[[
        "protein_accession_clean",
        "gene_symbol",
        "protein_name_final"
    ]].drop_duplicates()

    template["process_label"] = ""
    template["process_group"] = ""
    template["annotation_source"] = ""
    template["evidence_note"] = ""

    template.to_csv(OUTFILE, index=False)
    print(f"Wrote: {OUTFILE}")
    print(template.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
