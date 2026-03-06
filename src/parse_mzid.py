
## 8) First parser script

#Put this in `src/parse_mzid.py`:

from pathlib import Path
import pandas as pd
from pyteomics import mzid


RAW_DIR = Path("data/mzmid")
OUT_DIR = Path("data/processed")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_one_mzid(mzid_path: Path) -> pd.DataFrame:
    rows = []

    with mzid.read(str(mzid_path), retrieve_refs=True, read_schema=False) as reader:
        for sir in reader:
            spectrum_id = sir.get("spectrumID")
            spectra_data_ref = sir.get("spectraData_ref")
            sii_list = sir.get("SpectrumIdentificationItem", [])

            for sii in sii_list:
                peptide_sequence = None
                peptide_obj = sii.get("Peptide")
                if isinstance(peptide_obj, dict):
                    peptide_sequence = peptide_obj.get("PeptideSequence")

                pe_refs = sii.get("PeptideEvidenceRef", [])
                if not isinstance(pe_refs, list):
                    pe_refs = [pe_refs]

                if not pe_refs:
                    rows.append({
                        "source_file": mzid_path.name,
                        "spectrum_id": spectrum_id,
                        "spectra_data_ref": spectra_data_ref,
                        "rank": sii.get("rank"),
                        "charge": sii.get("chargeState"),
                        "pass_threshold": sii.get("passThreshold"),
                        "calculated_mz": sii.get("calculatedMassToCharge"),
                        "experimental_mz": sii.get("experimentalMassToCharge"),
                        "peptide_sequence": peptide_sequence,
                        "protein_accession": None
                    })
                else:
                    for pe in pe_refs:
                        protein_accession = None
                        if isinstance(pe, dict):
                            dbseq = pe.get("DBSequence")
                            if isinstance(dbseq, dict):
                                protein_accession = dbseq.get("accession")

                        rows.append({
                            "source_file": mzid_path.name,
                            "spectrum_id": spectrum_id,
                            "spectra_data_ref": spectra_data_ref,
                            "rank": sii.get("rank"),
                            "charge": sii.get("chargeState"),
                            "pass_threshold": sii.get("passThreshold"),
                            "calculated_mz": sii.get("calculatedMassToCharge"),
                            "experimental_mz": sii.get("experimentalMassToCharge"),
                            "peptide_sequence": peptide_sequence,
                            "protein_accession": protein_accession
                        })

    return pd.DataFrame(rows)


def main():
    all_dfs = []

    for mzid_file in sorted(RAW_DIR.glob("*.mzid")):
        print(f"Parsing {mzid_file.name}")
        df = parse_one_mzid(mzid_file)
        all_dfs.append(df)

    if not all_dfs:
        raise FileNotFoundError(f"No mzid files found in {RAW_DIR}")

    psm_df = pd.concat(all_dfs, ignore_index=True)
    psm_df.to_csv(OUT_DIR / "psm_table.csv", index=False)

    protein_df = (
        psm_df
        .dropna(subset=["protein_accession"])
        .groupby(["source_file", "protein_accession"], as_index=False)
        .agg(
            n_psms=("spectrum_id", "count"),
            n_peptides=("peptide_sequence", "nunique"),
            best_rank=("rank", "min"),
            any_pass_threshold=("pass_threshold", lambda x: any(str(v).lower() == "true" for v in x))
        )
        .sort_values(["source_file", "n_peptides", "n_psms"], ascending=[True, False, False])
    )

    protein_df.to_csv(OUT_DIR / "protein_evidence_summary.csv", index=False)

    print("Wrote:")
    print(OUT_DIR / "psm_table.csv")
    print(OUT_DIR / "protein_evidence_summary.csv")


if __name__ == "__main__":
    main()
