#!/usr/bin/env python3

from pathlib import Path
import re
import time
import requests
import pandas as pd

INFILE = Path("data/processed/protein_evidence_summary.csv")
OUTFILE = Path("data/processed/protein_evidence_annotated.csv")
MAPFILE = Path("data/processed/uniprot_annotation_map.csv")

# Optional local fallback FASTA
LOCAL_FASTA = Path("data/metadata/uniprot*.fasta")

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 25
SLEEP_SEC = 0.5
REQUEST_TIMEOUT = 60


def chunked(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def clean_accession(acc: str) -> str:
    if pd.isna(acc):
        return None
    acc = str(acc).strip()
    if acc.startswith("CON__"):
        acc = acc.replace("CON__", "", 1)
    return acc or None


def is_contaminant(acc: str) -> bool:
    if pd.isna(acc):
        return False
    return str(acc).startswith("CON__")


def parse_fasta_header(header: str):
    """
    Parse common UniProt FASTA header forms, e.g.
    >sp|P31327|CFL1_HUMAN Cofilin-1 OS=Homo sapiens OX=9606 GN=CFL1 PE=1 SV=3
    >tr|A0A0...|... Protein name OS=... GN=...
    """
    header = header.strip()
    if header.startswith(">"):
        header = header[1:]

    parts = header.split("|")
    accession = None
    protein_name = None
    gene_symbol = None
    organism_name = None
    uniprot_id = None

    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        accession = parts[1].strip()
        rest = parts[2]

        # rest begins like "CFL1_HUMAN Cofilin-1 OS=Homo sapiens ..."
        rest_parts = rest.split(" ", 1)
        if len(rest_parts) == 2:
            uniprot_id = rest_parts[0].strip()
            desc = rest_parts[1]
        else:
            desc = rest

        protein_name_match = re.split(r"\sOS=", desc, maxsplit=1)
        protein_name = protein_name_match[0].strip() if protein_name_match else None

        gn_match = re.search(r"\bGN=([^\s]+)", desc)
        if gn_match:
            gene_symbol = gn_match.group(1)

        os_match = re.search(r"\bOS=(.+?)\sOX=", desc)
        if os_match:
            organism_name = os_match.group(1).strip()

    return {
        "protein_accession_clean": accession,
        "uniprot_id": uniprot_id,
        "gene_symbol": gene_symbol,
        "protein_name_uniprot": protein_name,
        "organism_name": organism_name,
        "annotation_source": "local_fasta"
    }


def load_local_fasta_annotations(fasta_path: Path) -> pd.DataFrame:
    if not fasta_path.exists():
        return pd.DataFrame(columns=[
            "protein_accession_clean",
            "uniprot_id",
            "gene_symbol",
            "protein_name_uniprot",
            "organism_name",
            "annotation_source"
        ])

    rows = []
    with open(fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                parsed = parse_fasta_header(line.rstrip())
                if parsed["protein_accession_clean"]:
                    rows.append(parsed)

    if not rows:
        return pd.DataFrame(columns=[
            "protein_accession_clean",
            "uniprot_id",
            "gene_symbol",
            "protein_name_uniprot",
            "organism_name",
            "annotation_source"
        ])

    return pd.DataFrame(rows).drop_duplicates(subset=["protein_accession_clean"])


def extract_uniprot_result(result: dict):
    accession = result.get("primaryAccession")
    uniprot_id = result.get("uniProtkbId")

    protein_name = None
    protein_desc = result.get("proteinDescription", {})
    recommended = protein_desc.get("recommendedName", {})
    if isinstance(recommended, dict):
        full_name = recommended.get("fullName", {})
        if isinstance(full_name, dict):
            protein_name = full_name.get("value")

    if protein_name is None:
        submission = protein_desc.get("submissionNames", [])
        if submission and isinstance(submission, list):
            first = submission[0]
            if isinstance(first, dict):
                full_name = first.get("fullName", {})
                if isinstance(full_name, dict):
                    protein_name = full_name.get("value")

    gene_symbol = None
    genes = result.get("genes", [])
    if genes and isinstance(genes, list):
        first_gene = genes[0]
        if isinstance(first_gene, dict):
            gene_name_obj = first_gene.get("geneName", {})
            if isinstance(gene_name_obj, dict):
                gene_symbol = gene_name_obj.get("value")

    organism_name = result.get("organism", {}).get("scientificName")

    return {
        "protein_accession_clean": accession,
        "uniprot_id": uniprot_id,
        "gene_symbol": gene_symbol,
        "protein_name_uniprot": protein_name,
        "organism_name": organism_name,
        "annotation_source": "uniprot_api"
    }


def fetch_batch(batch):
    """
    Query UniProt for a batch of accessions.
    Uses small OR queries; if this fails, caller can fall back to single queries.
    """
    query = " OR ".join(batch)

    params = {
        "query": query,
        "format": "json",
        "fields": "accession,id,protein_name,gene_primary,organism_name",
        "size": max(len(batch), 25)
    }

    r = requests.get(UNIPROT_API, params=params, timeout=REQUEST_TIMEOUT)

    if not r.ok:
        raise requests.HTTPError(
            f"Status {r.status_code}. Response: {r.text[:500]}",
            response=r
        )

    data = r.json()
    rows = [extract_uniprot_result(result) for result in data.get("results", [])]
    return rows


def fetch_one(accession):
    """
    Fallback: query UniProt with a single accession.
    """
    params = {
        "query": accession,
        "format": "json",
        "fields": "accession,id,protein_name,gene_primary,organism_name",
        "size": 5
    }

    r = requests.get(UNIPROT_API, params=params, timeout=REQUEST_TIMEOUT)
    if not r.ok:
        return None

    data = r.json()
    results = data.get("results", [])
    if not results:
        return None

    # Prefer exact accession match
    for result in results:
        if result.get("primaryAccession") == accession:
            return extract_uniprot_result(result)

    return extract_uniprot_result(results[0])


def fetch_uniprot_annotations(accessions):
    all_rows = []

    for i, batch in enumerate(chunked(accessions, BATCH_SIZE), start=1):
        print(f"Batch {i}: querying {len(batch)} accessions")
        try:
            rows = fetch_batch(batch)
            all_rows.extend(rows)
        except Exception as e:
            print(f"  Batch failed: {e}")
            print("  Falling back to one accession at a time...")
            for acc in batch:
                row = fetch_one(acc)
                if row is not None:
                    all_rows.append(row)
                else:
                    all_rows.append({
                        "protein_accession_clean": acc,
                        "uniprot_id": None,
                        "gene_symbol": None,
                        "protein_name_uniprot": None,
                        "organism_name": None,
                        "annotation_source": "uniprot_api_failed"
                    })
                time.sleep(SLEEP_SEC)

        time.sleep(SLEEP_SEC)

    if not all_rows:
        return pd.DataFrame(columns=[
            "protein_accession_clean",
            "uniprot_id",
            "gene_symbol",
            "protein_name_uniprot",
            "organism_name",
            "annotation_source"
        ])

    return pd.DataFrame(all_rows).drop_duplicates(subset=["protein_accession_clean"])


def main():
    if not INFILE.exists():
        raise FileNotFoundError(f"Missing input file: {INFILE}")

    df = pd.read_csv(INFILE)

    df["is_contaminant"] = df["protein_accession"].apply(is_contaminant)
    df["protein_accession_clean"] = df["protein_accession"].apply(clean_accession)

    # Only annotate non-contaminants with non-empty accessions
    df_clean = df.loc[
        (~df["is_contaminant"]) &
        (df["protein_accession_clean"].notna())
    ].copy()

    accessions = sorted(df_clean["protein_accession_clean"].dropna().unique().tolist())
    print(f"Unique non-contaminant accessions to annotate: {len(accessions):,}")

    # Load optional local FASTA map first
    local_ann = load_local_fasta_annotations(LOCAL_FASTA)
    if not local_ann.empty:
        print(f"Loaded local FASTA annotations: {len(local_ann):,}")

    local_accs = set(local_ann["protein_accession_clean"].dropna().tolist()) if not local_ann.empty else set()
    missing_for_api = [acc for acc in accessions if acc not in local_accs]

    print(f"Accessions still needing UniProt API lookup: {len(missing_for_api):,}")

    api_ann = fetch_uniprot_annotations(missing_for_api) if missing_for_api else pd.DataFrame(columns=[
        "protein_accession_clean",
        "uniprot_id",
        "gene_symbol",
        "protein_name_uniprot",
        "organism_name",
        "annotation_source"
    ])

    ann = pd.concat([local_ann, api_ann], ignore_index=True) if not local_ann.empty else api_ann
    ann = ann.drop_duplicates(subset=["protein_accession_clean"], keep="first")
    ann.to_csv(MAPFILE, index=False)

    annotated = df.merge(ann, on="protein_accession_clean", how="left")
    annotated["protein_name_final"] = annotated["protein_name_uniprot"].fillna(annotated.get("protein_name"))

    preferred_cols = [
        "source_file",
        "protein_accession",
        "protein_accession_clean",
        "uniprot_id",
        "gene_symbol",
        "protein_name_final",
        "organism_name",
        "annotation_source",
        "n_psms",
        "n_peptides",
        "best_rank",
        "any_pass_threshold",
        "any_decoy",
        "is_contaminant"
    ]
    existing_cols = [c for c in preferred_cols if c in annotated.columns]
    other_cols = [c for c in annotated.columns if c not in existing_cols]
    annotated = annotated[existing_cols + other_cols]

    annotated.to_csv(OUTFILE, index=False)

    print(f"\nWrote: {MAPFILE}")
    print(f"Wrote: {OUTFILE}")

    print("\nAnnotation coverage:")
    print("Rows:", len(annotated))
    print("Annotated gene_symbol:", int(annotated["gene_symbol"].notna().sum()))
    print("Annotated protein_name_final:", int(annotated["protein_name_final"].notna().sum()))
    print("Contaminants flagged:", int(annotated["is_contaminant"].sum()))

    print("\nTop annotated rows:")
    show_cols = [
        "source_file", "protein_accession", "protein_accession_clean",
        "gene_symbol", "protein_name_final", "n_peptides", "is_contaminant",
        "annotation_source"
    ]
    show_cols = [c for c in show_cols if c in annotated.columns]
    print(annotated[show_cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
