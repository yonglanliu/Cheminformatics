# Copyright (c) 2026 Yonglan Liu
# Licensed under the MIT License.

import argparse
import csv
import io
import logging
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


# -----------------------------
# Logging
# -----------------------------
def setup_logger(log_path: str, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger("pubchem_pipeline")
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
# Config
# -----------------------------
@dataclass
class PipelineConfig:
    uniprot_id: str
    work_dir: Path = Path("./temp_pubchem")
    out_final_csv: Optional[Path] = None
    out_bio_csv: Optional[Path] = None
    log_path: str = "pubchem_pipeline.log"
    pause_sec: float = 0.25
    batch_size: int = 200
    properties: str = "CanonicalSMILES,IsomericSMILES,InChI,InChIKey"
    curl_retries: int = 4
    curl_timeout_sec: int = 180


# -----------------------------
# Helpers
# -----------------------------
def pick_col(df: pd.DataFrame, candidates: List[str]) -> str:
    cols = list(df.columns)
    lower_map = {c.lower(): c for c in cols}
    for c in candidates:
        if c in cols:
            return c
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    raise KeyError(f"Could not find any of columns: {candidates}. Example columns: {cols[:40]}")


def chunked(xs: List[int], n: int):
    for i in range(0, len(xs), n):
        yield xs[i : i + n]


def unit_to_molar_factor(u: Any) -> Optional[float]:
    """
    Convert an activity unit to a factor that turns Value * factor = M.
    Supports common units. Returns None if unknown.
    """
    if pd.isna(u):
        return None
    s = str(u).strip().replace("μ", "u").replace("µ", "u").lower().replace(" ", "")
    if s in ("m", "molar", "mol/l", "mol/litre"):
        return 1.0
    return {"pm": 1e-12, "nm": 1e-9, "um": 1e-6, "mm": 1e-3}.get(s)


def compute_pic50_from_value_and_unit(value: pd.Series, unit: pd.Series) -> pd.Series:
    """
    Unit-aware pIC50 computation:
      IC50_M = value * factor(unit)
      pIC50 = -log10(IC50_M)
    """
    factors = unit.apply(unit_to_molar_factor)
    ic50_m = value * factors
    pic50 = -np.log10(ic50_m)
    return pic50


# -----------------------------
# Curl wrapper with retry/backoff
# -----------------------------
def run_curl(
    url: str,
    *,
    out_path: Optional[Path],
    logger: logging.Logger,
    pause_sec: float,
    retries: int,
    timeout_sec: int,
) -> str:
    """
    If out_path is provided, downloads to file.
    Otherwise returns stdout.
    Retries with exponential backoff on failures or empty output.
    """
    cmd = ["curl", "-sSL", "--max-time", str(timeout_sec), url]

    for attempt in range(1, retries + 1):
        try:
            if out_path is not None:
                out_path.parent.mkdir(parents=True, exist_ok=True)
                full_cmd = cmd + ["-o", str(out_path)]
                logger.info(f"curl (file) attempt {attempt}/{retries}: {url}")
                subprocess.run(full_cmd, check=True)
                time.sleep(pause_sec)

                # sanity check file size
                if out_path.exists() and out_path.stat().st_size > 0:
                    return ""
                raise RuntimeError(f"Downloaded file is empty: {out_path}")

            else:
                logger.info(f"curl (stdout) attempt {attempt}/{retries}: {url}")
                res = subprocess.run(cmd, check=True, capture_output=True, text=True)
                time.sleep(pause_sec)

                if res.stdout and res.stdout.strip():
                    return res.stdout
                raise RuntimeError("Empty stdout from curl")

        except Exception as e:
            wait = min(8.0, 0.5 * (2 ** (attempt - 1)))
            logger.warning(f"curl failed: {e}. Backing off {wait:.1f}s")
            time.sleep(wait)

    raise RuntimeError(f"curl failed after {retries} attempts: {url}")


# -----------------------------
# PubChem fetchers
# -----------------------------
def fetch_target_concise_csv(cfg: PipelineConfig, logger: logging.Logger) -> Path:
    out_csv = cfg.work_dir / f"{cfg.uniprot_id}_target_concise.csv"
    url = f"{BASE}/assay/target/accession/{cfg.uniprot_id}/concise/CSV"
    run_curl(
        url,
        out_path=out_csv,
        logger=logger,
        pause_sec=cfg.pause_sec,
        retries=cfg.curl_retries,
        timeout_sec=cfg.curl_timeout_sec,
    )
    return out_csv


def fetch_structure_on_cid(cfg: PipelineConfig, cids: Iterable[int], logger: logging.Logger) -> pd.DataFrame:
    cid_list = sorted({int(c) for c in cids if pd.notna(c)})
    if not cid_list:
        return pd.DataFrame(columns=["CID"] + cfg.properties.split(","))

    frames = []
    for i, batch in enumerate(chunked(cid_list, cfg.batch_size), start=1):
        cid_str = ",".join(map(str, batch))
        url = f"{BASE}/compound/cid/{cid_str}/property/{cfg.properties}/CSV"

        logger.info(f"Fetching structures batch {i} ({len(batch)} CIDs)")
        text = run_curl(
            url,
            out_path=None,
            logger=logger,
            pause_sec=cfg.pause_sec,
            retries=cfg.curl_retries,
            timeout_sec=cfg.curl_timeout_sec,
        )

        df_props = pd.read_csv(io.StringIO(text))
        df_props["CID"] = pd.to_numeric(df_props["CID"], errors="coerce").astype("Int64")
        frames.append(df_props)

    out = pd.concat(frames, ignore_index=True).dropna(subset=["CID"])
    out["CID"] = out["CID"].astype(int)
    out = out.drop_duplicates(subset=["CID"])
    logger.info(f"Structures fetched: {len(out)} unique CIDs")
    return out


# -----------------------------
# Pipeline
# -----------------------------
def run_pipeline(cfg: PipelineConfig) -> pd.DataFrame:
    logger = setup_logger(cfg.log_path)
    logger.info("=== PubChem pipeline start ===")
    logger.info(f"UniProt: {cfg.uniprot_id}")
    logger.info(f"Work dir: {cfg.work_dir.resolve()}")

    cfg.work_dir.mkdir(parents=True, exist_ok=True)

    # 1) Download concise CSV
    logger.info("[Step 1] Download target concise CSV")
    concise_path = fetch_target_concise_csv(cfg, logger)
    logger.info(f"Downloaded: {concise_path} ({concise_path.stat().st_size} bytes)")

    if cfg.out_bio_csv:
        cfg.out_bio_csv.parent.mkdir(parents=True, exist_ok=True)
        cfg.out_bio_csv.write_bytes(concise_path.read_bytes())
        logger.info(f"Copied raw concise CSV to: {cfg.out_bio_csv}")

    # 2) Read
    logger.info("[Step 2] Read CSV into DataFrame")
    df = pd.read_csv(concise_path)
    logger.info(f"Rows: {len(df)} | Cols: {len(df.columns)}")

    # 3) Resolve columns
    logger.info("[Step 3] Resolve key columns")
    col_target = pick_col(df, ["Target Accession", "Protein Accession", "ProteinAccession", "Target.ProteinAccession"])
    col_a_name = pick_col(df, ["Activity Name", "Activity.Name", "Activity_Name", "Activity"])
    col_cid = pick_col(df, ["CID"])

    # value column candidates (PubChem varies)
    value_candidates = ["Activity Value [nM]", "Activity Value [uM]", "Activity Value", "Activity.Value", "Activity_Value", "Value"]
    unit_candidates = ["Activity Unit", "Activity.Unit", "Activity_Unit", "Unit"]

    col_value = None
    for c in value_candidates:
        if c in df.columns:
            col_value = c
            break
        # case-insensitive
        if c.lower() in {x.lower() for x in df.columns}:
            col_value = pick_col(df, [c])
            break
    if col_value is None:
        raise KeyError(f"Cannot find activity value column. Tried: {value_candidates}")

    try:
        col_unit = pick_col(df, unit_candidates)
    except KeyError:
        col_unit = None

    logger.info(f"Using columns: target='{col_target}', activity='{col_a_name}', CID='{col_cid}', value='{col_value}', unit='{col_unit}'")

    # 4) Filter strict target
    logger.info("[Step 4] Filter strict target (remove multi-target noise)")
    before = len(df)
    df = df[df[col_target].astype(str).str.strip() == cfg.uniprot_id].copy()
    logger.info(f"Rows: {before} -> {len(df)}")

    # 5) IC50 only
    logger.info("[Step 5] Filter IC50 only")
    before = len(df)
    df = df[df[col_a_name].astype(str).str.strip().str.upper() == "IC50"].copy()
    logger.info(f"Rows: {before} -> {len(df)}")

    # 6) Clean CID/value
    logger.info("[Step 6] Clean CID and numeric values")
    before = len(df)
    df[col_cid] = pd.to_numeric(df[col_cid], errors="coerce")
    df[col_value] = pd.to_numeric(df[col_value], errors="coerce")
    df = df.dropna(subset=[col_cid, col_value]).copy()
    df[col_cid] = df[col_cid].astype(int)
    logger.info(f"Rows: {before} -> {len(df)}")

    # 7) Compute pIC50 (nM vs uM vs unit-aware)
    logger.info("[Step 7] Compute pIC50")
    if "[nM]" in col_value:
        # nM -> pIC50 = 9 - log10(nM)
        df = df[df[col_value] > 0].copy()
        df["pIC50"] = 9.0 - np.log10(df[col_value].astype(float))
        logger.info("Computed pIC50 assuming values are in nM.")
    elif "[uM]" in col_value:
        # uM -> pIC50 = 6 - log10(uM)
        df = df[df[col_value] > 0].copy()
        df["pIC50"] = 6.0 - np.log10(df[col_value].astype(float))
        logger.info("Computed pIC50 assuming values are in uM.")
    else:
        if col_unit is None:
            raise KeyError("Value column has no explicit unit tag and no unit column found.")
        df = df[df[col_value] > 0].copy()
        df["pIC50"] = compute_pic50_from_value_and_unit(df[col_value].astype(float), df[col_unit])
        before = len(df)
        df = df.dropna(subset=["pIC50"]).copy()
        logger.info(f"Computed unit-aware pIC50. Dropped rows with unknown units: {before} -> {len(df)}")

    df["pIC50"] = df["pIC50"].round(2)
    logger.info(f"Rows after pIC50: {len(df)}")

    # 8) Deduplicate per CID
    logger.info("[Step 8] Deduplicate per CID (median pIC50)")
    df_bio = df.groupby(col_cid, as_index=False)["pIC50"].median()
    
    df_bio.to_csv(cfg.out_bio_csv, index=False)
    logger.info(f"Unique CIDs: {len(df_bio)}")

    # 9) Fetch structures
    logger.info("[Step 9] Fetch SMILES for CIDs")
    df_struct = fetch_structure_on_cid(cfg, df_bio["CID"].tolist(), logger)

    smiles_col = "IsomericSMILES" if "IsomericSMILES" in df_struct.columns else "CanonicalSMILES"
    df_struct = df_struct.rename(columns={smiles_col: "SMILES"})
    df_struct = df_struct[["CID", "SMILES"]].dropna(subset=["SMILES"])
    logger.info(f"SMILES available for: {len(df_struct)} CIDs")

    # 10) Merge
    logger.info("[Step 10] Merge bioactivity + structure")
    df_out = df_bio.merge(df_struct, on="CID", how="inner")
    df_out = df_out[["CID", "SMILES", "pIC50"]].sort_values("CID").reset_index(drop=True)
    logger.info(f"Final rows: {len(df_out)}")

    # 11) Write
    if cfg.out_final_csv:
        cfg.out_final_csv.parent.mkdir(parents=True, exist_ok=True)
        df_out.to_csv(cfg.out_final_csv, index=False)
        logger.info(f"Wrote final CSV: {cfg.out_final_csv}")

    logger.info("=== PubChem pipeline done ===")
    return df_out


# -----------------------------
# CLI
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Fetch PubChem IC50 for UniProt and output SMILES + median pIC50.")
    parser.add_argument("--uniprot", required=True, help="UniProt accession (e.g., P00533 for EGFR)")
    parser.add_argument("--out_final_csv", default=None, help="Output CSV path (CID, SMILES, pIC50)")
    parser.add_argument("--out_bio_csv", default=None, help="Optional: save the raw target concise CSV here")
    parser.add_argument("--work_dir", default="./temp_pubchem", help="Working directory")
    parser.add_argument("--log_path", default="pubchem_pipeline.log", help="Log file path")
    parser.add_argument("--pause_sec", type=float, default=0.25, help="Polite delay between PubChem requests")
    parser.add_argument("--batch_size", type=int, default=200, help="CID batch size per structure request")

    args = parser.parse_args()

    cfg = PipelineConfig(
        uniprot_id=args.uniprot,
        work_dir=Path(args.work_dir),
        out_final_csv=Path(args.out_final_csv) if args.out_final_csv else None,
        out_bio_csv=Path(args.out_bio_csv) if args.out_bio_csv else None,
        log_path=args.log_path,
        pause_sec=args.pause_sec,
        batch_size=args.batch_size,
    )

    run_pipeline(cfg)


if __name__ == "__main__":
    main()
