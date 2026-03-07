#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

INFILE = Path("data/metadata/protein_process_annotations_template.csv")
OUTFILE = Path("data/metadata/protein_process_annotations_bootstrap.csv")

RULES = [
    ("receptor", "receptor-mediated endocytosis", "uptake"),
    ("lipoprotein receptor", "lipoprotein particle binding", "uptake"),
    ("transferrin receptor", "iron uptake", "iron_homeostasis"),
    ("hemopexin", "heme scavenging", "heme_handling"),
    ("haptoglobin", "heme/hemoglobin scavenging", "heme_handling"),
    ("ceruloplasmin", "iron homeostasis", "iron_homeostasis"),
    ("coagulation", "coagulation", "coagulation"),
    ("complement", "complement activation", "complement"),
    ("ficolin", "lectin pathway of complement activation", "complement"),
    ("apolipoprotein", "lipoprotein metabolism", "lipoprotein_biology"),
    ("serpin", "protease inhibition", "protease_regulation"),
    ("macroglobulin", "protease inhibition", "protease_regulation"),
    ("immunoglobulin", "immunoglobulin complex", "immune_background"),
]


MANUAL_OVERRIDES = {
    "LRP1": {
        "process_label": "receptor-mediated endocytosis;lipoprotein particle binding",
        "process_group": "uptake;uptake",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "Top HPX-specific candidate"
    },
    "TFRC": {
        "process_label": "iron uptake;receptor-mediated endocytosis",
        "process_group": "iron_homeostasis;uptake",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "Shared HPX/Tf uptake-related candidate"
    },
    "CP": {
        "process_label": "iron homeostasis;oxidation-reduction process",
        "process_group": "iron_homeostasis;iron_homeostasis",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "HPX-associated extracellular protein"
    },
    "HP": {
        "process_label": "heme/hemoglobin scavenging",
        "process_group": "heme_handling",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "HPX-associated extracellular protein"
    },
    "HPR": {
        "process_label": "heme/hemoglobin scavenging",
        "process_group": "heme_handling",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "HPX-associated extracellular protein"
    },
    "F12": {
        "process_label": "coagulation;proteolytic cascade",
        "process_group": "coagulation;coagulation",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "HPX-associated extracellular protein"
    },
    "SERPINA1": {
        "process_label": "protease inhibition",
        "process_group": "protease_regulation",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "HPX-associated extracellular protein"
    },
    "A2M": {
        "process_label": "protease inhibition;extracellular protein binding",
        "process_group": "protease_regulation;protease_regulation",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "HPX-associated extracellular protein"
    },
    "FCN3": {
        "process_label": "innate immune response;lectin pathway of complement activation",
        "process_group": "complement;complement",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "HPX-associated extracellular protein"
    },
    "H6PD": {
        "process_label": "glucose metabolism",
        "process_group": "metabolism",
        "annotation_source": "bootstrap_manual_override",
        "evidence_note": "HPX-associated protein"
    }
}


def infer_from_name(protein_name):
    if pd.isna(protein_name):
        return "", ""

    name = str(protein_name).lower()
    labels = []
    groups = []

    for keyword, label, group in RULES:
        if keyword in name:
            labels.append(label)
            groups.append(group)

    # de-duplicate while preserving order
    seen = set()
    dedup_labels = []
    dedup_groups = []
    for label, group in zip(labels, groups):
        key = (label, group)
        if key not in seen:
            seen.add(key)
            dedup_labels.append(label)
            dedup_groups.append(group)

    return ";".join(dedup_labels), ";".join(dedup_groups)


def main():
    if not INFILE.exists():
        raise FileNotFoundError(f"Missing input file: {INFILE}")

    df = pd.read_csv(INFILE)

    out = df.copy()

    # Fill from keyword rules
    inferred = out["protein_name_final"].apply(infer_from_name)
    out[["process_label", "process_group"]] = pd.DataFrame(inferred.tolist(), index=out.index)

    out["annotation_source"] = "bootstrap_name_keyword"
    out["evidence_note"] = "Auto-filled from protein_name keyword rules"

    # Apply manual overrides for key proteins
    for gene, vals in MANUAL_OVERRIDES.items():
        mask = out["gene_symbol"] == gene
        if mask.any():
            for col, val in vals.items():
                out.loc[mask, col] = val

    out.to_csv(OUTFILE, index=False)

    print(f"Wrote: {OUTFILE}")
    print("\nPreview:")
    print(out.head(20).to_string(index=False))

    print("\nRows with non-empty process_label:")
    print((out["process_label"].fillna("").str.strip() != "").sum())


if __name__ == "__main__":
    main()
