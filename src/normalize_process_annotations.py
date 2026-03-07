#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

INFILE = Path("data/metadata/protein_process_annotations_bootstrap.csv")
OUTFILE = Path("data/metadata/protein_process_annotations.csv")


def split_multi_value(val):
    if pd.isna(val):
        return []
    val = str(val).strip()
    if not val:
        return []
    return [x.strip() for x in val.split(";") if x.strip()]


def main():
    if not INFILE.exists():
        raise FileNotFoundError(f"Missing input file: {INFILE}")

    df = pd.read_csv(INFILE)

    rows = []
    for _, row in df.iterrows():
        process_labels = split_multi_value(row.get("process_label", ""))
        process_groups = split_multi_value(row.get("process_group", ""))

        if not process_labels:
            continue

        if len(process_groups) == 1 and len(process_labels) > 1:
            process_groups = process_groups * len(process_labels)

        if len(process_groups) < len(process_labels):
            process_groups = process_groups + [""] * (len(process_labels) - len(process_groups))

        for label, group in zip(process_labels, process_groups):
            rows.append({
                "protein_accession_clean": row.get("protein_accession_clean"),
                "gene_symbol": row.get("gene_symbol"),
                "protein_name_final": row.get("protein_name_final"),
                "process_label": label,
                "process_group": group,
                "annotation_source": row.get("annotation_source", ""),
                "evidence_note": row.get("evidence_note", "")
            })

    out = pd.DataFrame(rows).drop_duplicates()
    out.to_csv(OUTFILE, index=False)

    print(f"Wrote: {OUTFILE}")
    print("\nPreview:")
    print(out.head(25).to_string(index=False))


if __name__ == "__main__":
    main()
