"""
Created by: Yonglan Liu
Date: 2026-02-13
Purpose: Combine ChEMBL, PubChem, BindingDB datasets using structure keys (InChIKey)
         and track progress with a logger.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi


# -----------------------------
# Logging
# -----------------------------
def setup_logger(log_path: str, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger("combine_datasets")
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    # Console
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    # File
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# -----------------------------
# Chemistry helpers
# -----------------------------
def clean_smiles(s: object) -> Optional[str]:
    """Strip CXSMILES extensions like 'SMILES |...|' and whitespace."""
    if pd.isna(s):
        return None
    s = str(s).strip()
    if not s:
        return None
    if "|" in s:
        s = s.split("|", 1)[0].strip()
    return s or None


def smiles_to_inchikey(smiles: Optional[str]) -> Optional[str]:
    """
    Convert SMILES -> InChIKey.
    NOTE: Use MolToInchiKey (not MolToInchi) for merging keys.
    """
    if smiles is None:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


# -----------------------------
# Data processing
# -----------------------------
def _validate_cols(df: pd.DataFrame, needed: Tuple[str, ...], name: str, logger: logging.Logger) -> None:
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise KeyError(f"{name}: missing columns {missing}. Found columns: {list(df.columns)}")
    logger.info(f"{name}: columns OK ({', '.join(needed)})")


def _prep_df(path: str, source: str, logger: logging.Logger) -> pd.DataFrame:
    logger.info(f"Reading {source}: {path}")
    df = pd.read_csv(path)
    logger.info(f"{source}: loaded rows={len(df):,} cols={len(df.columns)}")

    _validate_cols(df, ("SMILES", "pIC50"), source, logger)

    # Clean SMILES
    before = len(df)
    has_pipe = df["SMILES"].astype(str).str.contains(r"\|", regex=True, na=False).sum()
    df["SMILES"] = df["SMILES"].apply(clean_smiles)
    df = df.dropna(subset=["SMILES"]).copy()
    logger.info(f"{source}: CXSMILES rows (had '|...|')={has_pipe:,}; after SMILES clean {before:,}->{len(df):,}")

    # pIC50 numeric + sane
    before = len(df)
    df["pIC50"] = pd.to_numeric(df["pIC50"], errors="coerce")
    df = df.dropna(subset=["pIC50"]).copy()
    df = df[df["pIC50"] > 0].copy()
    logger.info(f"{source}: after pIC50 numeric/positive {before:,}->{len(df):,}")

    # InChIKey
    logger.info(f"{source}: computing InChIKey from SMILES (RDKit)")
    before = len(df)
    df["InChIKey"] = df["SMILES"].apply(smiles_to_inchikey)
    bad = df["InChIKey"].isna().sum()
    df = df.dropna(subset=["InChIKey"]).copy()
    logger.info(f"{source}: InChIKey missing={bad:,}; after drop {before:,}->{len(df):,}")

    df["source"] = source
    return df[["InChIKey", "SMILES", "pIC50", "source"]]


def run_combine(
    chembl_data_path: str,
    pubchem_data_path: str,
    bindingdb_data_path: str,
    *,
    out_combined_path: Optional[str] = None,
    out_final_path: Optional[str] = None,
    out_overlap_path: Optional[str] = None,
    log_path: str = "combine_pipeline.log",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    logger = setup_logger(log_path)
    logger.info("=== Combine pipeline start ===")

    df_chembl = _prep_df(chembl_data_path, "ChEMBL", logger)
    df_pubchem = _prep_df(pubchem_data_path, "PubChem", logger)
    df_bindingdb = _prep_df(bindingdb_data_path, "BindingDB", logger)

    logger.info(f"Data from ChEMBL:   {len(df_chembl):,}")
    logger.info(f"Data from PubChem:  {len(df_pubchem):,}")
    logger.info(f"Data from BindingDB:{len(df_bindingdb):,}")

    # Combine (all measurements, keep source)
    combined = pd.concat([df_chembl, df_pubchem, df_bindingdb], ignore_index=True)
    logger.info(f"Combined total rows (all measurements): {len(combined):,}")
    logger.info(f"Unique InChIKeys (all measurements): {combined['InChIKey'].nunique():,}")

    # Overlap table (counts per source per InChIKey)
    overlap = (
        combined.groupby(["InChIKey", "source"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )

    logger.info("Overlap columns present: " + ", ".join([c for c in overlap.columns if c != "InChIKey"]))
    logger.info("Top 5 overlap preview:\n" + overlap.head().to_string(index=False))

    # Final dataset: one row per molecule (median pIC50 across sources)
    final_dataset = (
        combined.groupby("InChIKey", as_index=False)
        .agg(
            SMILES=("SMILES", "first"),
            pIC50=("pIC50", "median"),
            n_measurements=("pIC50", "count"),
            n_sources=("source", "nunique"),
        )
        .sort_values(["n_sources", "n_measurements"], ascending=False)
        .reset_index(drop=True)
    )

    logger.info(f"Final dataset rows (unique molecules): {len(final_dataset):,}")
    logger.info(
        "Final pIC50 stats: "
        f"min={final_dataset['pIC50'].min():.2f}, "
        f"median={final_dataset['pIC50'].median():.2f}, "
        f"max={final_dataset['pIC50'].max():.2f}"
    )

    # Write outputs
    if out_combined_path:
        Path(out_combined_path).parent.mkdir(parents=True, exist_ok=True)
        combined.to_csv(out_combined_path, index=False)
        logger.info(f"Wrote combined (all measurements): {out_combined_path}")

    if out_final_path:
        Path(out_final_path).parent.mkdir(parents=True, exist_ok=True)
        final_dataset.to_csv(out_final_path, index=False)
        logger.info(f"Wrote final (unique molecules): {out_final_path}")

    if out_overlap_path:
        Path(out_overlap_path).parent.mkdir(parents=True, exist_ok=True)
        overlap.to_csv(out_overlap_path, index=False)
        logger.info(f"Wrote overlap table: {out_overlap_path}")

    logger.info("=== Combine pipeline done ===")
    return combined, final_dataset, overlap


def main():
    parser = argparse.ArgumentParser(description="Combine ChEMBL, PubChem, BindingDB datasets by InChIKey.")
    parser.add_argument("--chembl", required=True, help="Path to ChEMBL CSV (must contain SMILES,pIC50)")
    parser.add_argument("--pubchem", required=True, help="Path to PubChem CSV (must contain SMILES,pIC50)")
    parser.add_argument("--bindingdb", required=True, help="Path to BindingDB CSV (must contain SMILES,pIC50)")
    parser.add_argument("--out_combined", default=None, help="Output CSV for all measurements (with source)")
    parser.add_argument("--out_final", default="final_merged_unique_molecules.csv", help="Output CSV for unique molecules")
    parser.add_argument("--out_overlap", default="overlap_by_source.csv", help="Output CSV for overlap table")
    parser.add_argument("--log_path", default="combine_pipeline.log", help="Log file path")
    args = parser.parse_args()

    run_combine(
        chembl_data_path=args.chembl,
        pubchem_data_path=args.pubchem,
        bindingdb_data_path=args.bindingdb,
        out_combined_path=args.out_combined,
        out_final_path=args.out_final,
        out_overlap_path=args.out_overlap,
        log_path=args.log_path,
    )


if __name__ == "__main__":
    main()
