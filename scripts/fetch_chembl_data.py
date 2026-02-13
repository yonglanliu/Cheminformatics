import argparse
import logging
from pathlib import Path
from typing import Union, Optional

import numpy as np
import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm.auto import tqdm


# -----------------------------
# Logging setup
# -----------------------------
def setup_logger(log_path: str | Path) -> logging.Logger:
    """
    Create a logger that writes to both console and a log file.
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("Fetch Chembl Data Pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()  # avoid duplicated handlers if run in notebooks

    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # File handler
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    return logger


# -----------------------------
# ChEMBL API clients
# -----------------------------
target_api = new_client.target
compounds_api = new_client.molecule
bioactivity_api = new_client.activity
document_api = new_client.document


# -----------------------------
# Helper functions
# -----------------------------
def fetch_bioactivity_by_target(target_chembl_id: str, logger: logging.Logger):
    """
    Fetch bioactivity records for a target from ChEMBL.
    We request Binding + Functional assays, and only exact relations (i.e., '=')
    to avoid censored values like '>' or '<'.
    """
    logger.info(f"Step 1: Fetching bioactivity for target {target_chembl_id}")

    # NOTE: In ChEMBL, the relation field commonly appears as "standard_relation".
    # Some client versions expose it as "relation". We'll request both safely.
    bioactivity = bioactivity_api.filter(
        target_chembl_id=target_chembl_id,
        assay_type__in=["B", "F"],   # Binding + Functional
        standard_relation="=",       # exact values only
    ).only(
        "activity_id",
        "document_chembl_id",
        "assay_chembl_id",
        "assay_description",
        "assay_type",
        "molecule_chembl_id",
        "standard_type",
        "standard_units",
        "standard_value",
        "standard_relation",
        "target_chembl_id",
        "target_organism",
    )

    # This triggers evaluation and gives a real count
    n = len(bioactivity)
    logger.info(f"Fetched {n} bioactivity rows from ChEMBL.")
    return bioactivity


def get_doi_safe(chembl_doc_id: str, logger: logging.Logger) -> str | None:
    """
    Retrieve DOI for a ChEMBL document ID. Returns None if missing/unavailable.
    """
    if chembl_doc_id is None or (isinstance(chembl_doc_id, float) and np.isnan(chembl_doc_id)):
        return None
    try:
        doc = document_api.get(chembl_doc_id)
        return doc.get("doi")
    except Exception as e:
        logger.warning(f"Could not fetch DOI for document {chembl_doc_id}: {e}")
        return None


def extract_smiles(molecule_structures) -> str | None:
    """
    molecule_structures is typically a dict like {'canonical_smiles': '...'}.
    Return canonical_smiles if present.
    """
    if isinstance(molecule_structures, dict):
        return molecule_structures.get("canonical_smiles")
    return None


def compute_pIC50_from_nM(ic50_nM: pd.Series) -> pd.Series:
    """
    Convert IC50 in nM to pIC50.

    pIC50 = -log10(IC50 in molar)
          = -log10(IC50_nM * 1e-9)
          = 9 - log10(IC50_nM)
    """
    return 9 - np.log10(ic50_nM)


# -----------------------------
# Main pipeline
# -----------------------------
def run_pipeline(
    target_chembl_id: str,
    out_structure_plus_bioactivity: str | None,
    out_bio_path: str | None,
    log_path: str,
):
    logger = setup_logger(log_path)
    logger.info("=== ChEMBL IC50 -> pIC50 pipeline started ===")

    # Step 1: Fetch bioactivity
    bioactivity = fetch_bioactivity_by_target(target_chembl_id, logger)
    bio_df = pd.DataFrame.from_records(bioactivity)
    logger.info(f"Step 2: Bioactivity dataframe created. Shape = {bio_df.shape}")

    # Step 3: Attach DOI (optional but useful)
    logger.info("Step 3: Fetching DOI for each document_chembl_id (may take time)...")
    bio_df["doi"] = bio_df["document_chembl_id"].apply(lambda x: get_doi_safe(x, logger))

    # Step 4: Keep only IC50 rows (standard_type) and ensure correct units
    # NOTE: If you include 'Log IC50' you must treat it differently.
    # Here we keep only IC50 measured in nM so conversion is valid.
    logger.info("Step 4: Filtering to IC50 (nM) with exact relation '=' ...")

    bio_df["standard_type"] = bio_df["standard_type"].astype(str)
    bio_df["standard_units"] = bio_df["standard_units"].astype(str)

    bio_df = bio_df[
        (bio_df["standard_type"].str.upper() == "IC50") &
        (bio_df["standard_units"].str.lower() == "nm") &
        (bio_df["standard_relation"] == "=")
    ].copy()

    logger.info(f"After IC50+nM+relation filtering: {len(bio_df)} rows")

    # Step 5: Clean numeric standard_value
    logger.info("Step 5: Cleaning standard_value and computing pIC50...")

    n_missing = bio_df["standard_value"].isna().sum()
    logger.info(f"standard_value NaNs before cleaning: {n_missing}")

    bio_df = bio_df.dropna(subset=["standard_value"]).copy()
    bio_df["standard_value"] = pd.to_numeric(bio_df["standard_value"], errors="coerce")
    bio_df = bio_df.dropna(subset=["standard_value"]).copy()

    # Remove non-positive values (log10 undefined)
    n_nonpos = (bio_df["standard_value"] <= 0).sum()
    if n_nonpos > 0:
        logger.warning(f"Found {n_nonpos} non-positive IC50 values; removing them.")
        bio_df = bio_df[bio_df["standard_value"] > 0].copy()

    bio_df["pIC50"] = compute_pIC50_from_nM(bio_df["standard_value"]).round(2)

    logger.info(f"Computed pIC50. Remaining rows: {len(bio_df)}")

    # Optional: save cleaned bioactivity table
    if out_bio_path:
        out_bio_path = str(out_bio_path)
        Path(out_bio_path).parent.mkdir(parents=True, exist_ok=True)
        bio_df.to_csv(out_bio_path, index=False)
        logger.info(f"Saved cleaned bioactivity table to: {out_bio_path}")

    # Step 6: Aggregate duplicate measurements per molecule
    logger.info("Step 6: Aggregating duplicates per molecule_chembl_id using median pIC50...")

    bio_small = bio_df[["molecule_chembl_id", "pIC50"]].copy()
    bio_df_median = (
        bio_small
        .groupby("molecule_chembl_id", as_index=False)["pIC50"]
        .median()
    )
    logger.info(f"Before aggregation: {len(bio_small)} rows")
    logger.info(f"After  aggregation: {len(bio_df_median)} unique molecules")

    # Step 7: Download structures for those molecules
    logger.info("Step 7: Downloading molecule structures from ChEMBL...")

    mol_ids = bio_df_median["molecule_chembl_id"].tolist()
    compounds_provider = compounds_api.filter(
        molecule_chembl_id__in=mol_ids
    ).only("molecule_chembl_id", "molecule_structures")

    compounds = list(tqdm(compounds_provider, desc="Downloading compounds"))
    com_df = pd.DataFrame.from_records(compounds)
    logger.info(f"Downloaded {len(com_df)} compounds. Shape = {com_df.shape}")

    # Step 8: Extract SMILES and drop missing
    logger.info("Step 8: Extracting canonical SMILES and removing missing SMILES...")

    com_df["SMILES"] = com_df["molecule_structures"].apply(extract_smiles)
    missing_smiles = com_df["SMILES"].isna().sum()
    logger.info(f"{missing_smiles} compounds missing SMILES (will be removed).")

    com_df = com_df.dropna(subset=["SMILES"]).copy()
    logger.info(f"Compounds remaining after SMILES filter: {len(com_df)}")

    # Step 9: Merge median pIC50 with SMILES
    logger.info("Step 9: Merging median pIC50 with SMILES...")

    combine_df = bio_df_median.merge(
        com_df[["molecule_chembl_id", "SMILES"]],
        on="molecule_chembl_id",
        how="inner",
    )

    logger.info(f"Final merged dataset size: {len(combine_df)} rows")
    logger.info(f"Final columns: {list(combine_df.columns)}")

    # Save final dataset
    if out_structure_plus_bioactivity:
        out_structure_plus_bioactivity = str(out_structure_plus_bioactivity)
        Path(out_structure_plus_bioactivity).parent.mkdir(parents=True, exist_ok=True)
        combine_df.to_csv(out_structure_plus_bioactivity, index=False)
        logger.info(f"Saved final dataset to: {out_structure_plus_bioactivity}")

    logger.info("=== Pipeline finished successfully ===")


def main():
    parser = argparse.ArgumentParser(description="Fetch ChEMBL IC50 data and build SMILES+pIC50 dataset.")
    parser.add_argument("--target_chembl_id", required=True, help="Target ChEMBL ID, e.g. CHEMBLxxxx")
    parser.add_argument("--out_final_csv", default=None, help="Output CSV path for SMILES + median pIC50")
    parser.add_argument("--out_bio_csv", default=None, help="Output CSV path for cleaned bioactivity table")
    parser.add_argument("--log_path", default="chembl_pipeline.log", help="Log file path")

    args = parser.parse_args()

    run_pipeline(
        target_chembl_id=args.target_chembl_id,
        out_structure_plus_bioactivity=args.out_final_csv,
        out_bio_path=args.out_bio_csv,
        log_path=args.log_path,
    )


if __name__ == "__main__":
    main()
