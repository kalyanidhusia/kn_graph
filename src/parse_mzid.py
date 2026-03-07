#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import xml.etree.ElementTree as ET

RAW_DIR = Path("data/mzmid")
OUT_DIR = Path("data/processed")
OUT_DIR.mkdir(parents=True, exist_ok=True)

NS = {"mzid": "http://psidev.info/psi/pi/mzIdentML/1.1"}


def safe_bool(x):
    if pd.isna(x):
        return False
    return str(x).strip().lower() == "true"


def parse_one_mzid_xml(mzid_path: Path) -> pd.DataFrame:
    tree = ET.parse(mzid_path)
    root = tree.getroot()

    # Map DBSequence ids -> protein accession/name
    dbseq_map = {}
    for dbs in root.findall(".//mzid:DBSequence", NS):
        dbseq_id = dbs.attrib.get("id")
        accession = dbs.attrib.get("accession")
        name = dbs.attrib.get("name")
        search_db_ref = dbs.attrib.get("searchDatabase_ref")

        dbseq_map[dbseq_id] = {
            "accession": accession,
            "name": name,
            "search_db_ref": search_db_ref,
        }

    # Map Peptide ids -> peptide sequence
    peptide_map = {}
    for pep in root.findall(".//mzid:Peptide", NS):
        pep_id = pep.attrib.get("id")
        seq_elem = pep.find("mzid:PeptideSequence", NS)
        peptide_seq = seq_elem.text if seq_elem is not None else None
        peptide_map[pep_id] = peptide_seq

    # Map PeptideEvidence ids -> peptide_ref + DBSequence ref
    peptide_evidence_map = {}
    for pe in root.findall(".//mzid:PeptideEvidence", NS):
        pe_id = pe.attrib.get("id")
        peptide_evidence_map[pe_id] = {
            "peptide_ref": pe.attrib.get("peptide_ref"),
            "dbsequence_ref": pe.attrib.get("dBSequence_ref"),
            "start": pe.attrib.get("start"),
            "end": pe.attrib.get("end"),
            "pre": pe.attrib.get("pre"),
            "post": pe.attrib.get("post"),
            "is_decoy": pe.attrib.get("isDecoy"),
        }

    rows = []

    for sir in root.findall(".//mzid:SpectrumIdentificationResult", NS):
        spectrum_id = sir.attrib.get("spectrumID")
        spectra_data_ref = sir.attrib.get("spectraData_ref")
        result_id = sir.attrib.get("id")
        result_name = sir.attrib.get("name")

        spectrum_title = None
        for cv in sir.findall("mzid:cvParam", NS):
            if cv.attrib.get("name") == "spectrum title":
                spectrum_title = cv.attrib.get("value")

        for sii in sir.findall("mzid:SpectrumIdentificationItem", NS):
            sii_id = sii.attrib.get("id")
            rank = sii.attrib.get("rank")
            charge = sii.attrib.get("chargeState")
            pass_threshold = sii.attrib.get("passThreshold")
            calc_mz = sii.attrib.get("calculatedMassToCharge")
            exp_mz = sii.attrib.get("experimentalMassToCharge")
            peptide_ref = sii.attrib.get("peptide_ref")

            peptide_sequence = peptide_map.get(peptide_ref)

            # Capture all score-like cvParams
            score_dict = {}
            for cv in sii.findall("mzid:cvParam", NS):
                name = cv.attrib.get("name")
                value = cv.attrib.get("value")
                if name:
                    score_dict[name] = value

            pe_refs = sii.findall("mzid:PeptideEvidenceRef", NS)

            if not pe_refs:
                rows.append({
                    "source_file": mzid_path.name,
                    "sir_id": result_id,
                    "sir_name": result_name,
                    "spectrum_id": spectrum_id,
                    "spectrum_title": spectrum_title,
                    "spectra_data_ref": spectra_data_ref,
                    "sii_id": sii_id,
                    "rank": rank,
                    "charge": charge,
                    "pass_threshold": pass_threshold,
                    "calculated_mz": calc_mz,
                    "experimental_mz": exp_mz,
                    "peptide_ref": peptide_ref,
                    "peptide_sequence": peptide_sequence,
                    "peptide_evidence_ref": None,
                    "protein_accession": None,
                    "protein_name": None,
                    "protein_dbseq_ref": None,
                    "peptide_start": None,
                    "peptide_end": None,
                    "is_decoy": None,
                    **score_dict
                })
            else:
                for ref in pe_refs:
                    pe_id = ref.attrib.get("peptideEvidence_ref")
                    pe_info = peptide_evidence_map.get(pe_id, {})
                    dbseq_ref = pe_info.get("dbsequence_ref")
                    db_info = dbseq_map.get(dbseq_ref, {})

                    rows.append({
                        "source_file": mzid_path.name,
                        "sir_id": result_id,
                        "sir_name": result_name,
                        "spectrum_id": spectrum_id,
                        "spectrum_title": spectrum_title,
                        "spectra_data_ref": spectra_data_ref,
                        "sii_id": sii_id,
                        "rank": rank,
                        "charge": charge,
                        "pass_threshold": pass_threshold,
                        "calculated_mz": calc_mz,
                        "experimental_mz": exp_mz,
                        "peptide_ref": peptide_ref,
                        "peptide_sequence": peptide_sequence,
                        "peptide_evidence_ref": pe_id,
                        "protein_accession": db_info.get("accession"),
                        "protein_name": db_info.get("name"),
                        "protein_dbseq_ref": dbseq_ref,
                        "peptide_start": pe_info.get("start"),
                        "peptide_end": pe_info.get("end"),
                        "is_decoy": pe_info.get("is_decoy"),
                        **score_dict
                    })

    return pd.DataFrame(rows)


def build_protein_summary(psm_df: pd.DataFrame) -> pd.DataFrame:
    if psm_df.empty:
        return pd.DataFrame()

    # Optional protein name mapping if available
    protein_name_map = (
        psm_df[["protein_accession", "protein_name"]]
        .dropna(subset=["protein_accession"])
        .drop_duplicates()
    )

    protein_df = (
        psm_df
        .dropna(subset=["protein_accession"])
        .groupby(["source_file", "protein_accession"], as_index=False)
        .agg(
            n_psms=("spectrum_id", "count"),
            n_peptides=("peptide_sequence", pd.Series.nunique),
            best_rank=("rank", "min"),
            any_pass_threshold=("pass_threshold", lambda x: any(safe_bool(v) for v in x)),
            any_decoy=("is_decoy", lambda x: any(safe_bool(v) for v in x)),
        )
        .merge(protein_name_map, on="protein_accession", how="left")
        .sort_values(
            ["source_file", "n_peptides", "n_psms"],
            ascending=[True, False, False]
        )
    )

    return protein_df


def main():
    mzid_files = sorted(RAW_DIR.glob("*.mzid"))
    if not mzid_files:
        raise FileNotFoundError(f"No .mzid files found in {RAW_DIR}")

    print(f"Found {len(mzid_files)} mzID files in {RAW_DIR}")

    all_dfs = []
    for mzid_file in mzid_files:
        print(f"Parsing {mzid_file.name}")
        df = parse_one_mzid_xml(mzid_file)

        nonnull_pep = df["peptide_sequence"].notna().sum() if "peptide_sequence" in df.columns else 0
        nonnull_prot = df["protein_accession"].notna().sum() if "protein_accession" in df.columns else 0

        print(f"  rows={len(df):,} peptide_sequence_nonnull={nonnull_pep:,} protein_accession_nonnull={nonnull_prot:,}")
        all_dfs.append(df)

    psm_df = pd.concat(all_dfs, ignore_index=True)
    psm_out = OUT_DIR / "psm_table.csv"
    psm_df.to_csv(psm_out, index=False)

    protein_df = build_protein_summary(psm_df)
    protein_out = OUT_DIR / "protein_evidence_summary.csv"
    protein_df.to_csv(protein_out, index=False)

    print("\nFinished.")
    print(f"PSM table: {psm_out} ({len(psm_df):,} rows)")
    print(f"Protein summary: {protein_out} ({len(protein_df):,} rows)")

    if not protein_df.empty:
        print("\nTop protein summary rows:")
        print(protein_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
