"""
Created by: Yonglan Liu
Date: 2026-02-13
"""

import argparse
import json
import logging
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd


# -----------------------------
# Logging
# -----------------------------
def setup_logger(log_path: str, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger("bindingdb_pipeline")
    logger.setLevel(level)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    return logger


# -----------------------------
# Config
# -----------------------------
@dataclass
class BDBConfig:
    uniprot: str
    out_csv: Optional[Path] = None
    log_path: str = "bindingdb_pipeline.log"
    pause_sec: float = 0.25
    curl_retries: int = 4
    curl_timeout_sec: int = 180
    # You told me unit is uM:
    ic50_unit: str = "nM"   # "uM" or "nM" or "M"


# -----------------------------
# Curl helper (retry/backoff)
# -----------------------------
def curl_json(url: str, logger: logging.Logger, pause_sec: float, retries: int, timeout_sec: int) -> Dict[str, Any]:
    cmd = ["curl", "-sSL", "--max-time", str(timeout_sec), url]

    for attempt in range(1, retries + 1):
        try:
            logger.info(f"curl attempt {attempt}/{retries}: {url}")
            res = subprocess.run(cmd, check=True, capture_output=True, text=True)
            time.sleep(pause_sec)

            text = res.stdout.strip()
            if not text:
                raise RuntimeError("Empty response")

            return json.loads(text)

        except Exception as e:
            wait = min(8.0, 0.5 * (2 ** (attempt - 1)))
            logger.warning(f"curl/json failed: {e}. Backoff {wait:.1f}s")
            time.sleep(wait)

    raise RuntimeError(f"Failed after {retries} attempts: {url}")


# -----------------------------
# Core logic
# -----------------------------
def compute_pic50(series: pd.Series, unit: str) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    x = x.where(x > 0)

    u = unit.lower()
    if u == "um":
        return 6.0 - np.log10(x)
    if u == "nm":
        return 9.0 - np.log10(x)
    if u == "m":
        return -np.log10(x)

    raise ValueError(f"Unknown unit: {unit}")


def fetch_bindingdb_affinities(uniprot: str, logger: logging.Logger, cfg: BDBConfig) -> pd.DataFrame:
    # NOTE: your response key is "getLindsByUniprotResponse" (typo in API naming).
    url = f"https://bindingdb.org/rest/getLigandsByUniprot?uniprot={uniprot}&response=application/json"
    data = curl_json(url, logger, cfg.pause_sec, cfg.curl_retries, cfg.curl_timeout_sec)

    top_key = "getLindsByUniprotResponse"
    if top_key not in data:
        raise KeyError(f"Unexpected JSON keys: {list(data.keys())[:20]} (expected '{top_key}')")

    resp = data[top_key]

    hits = resp.get("bdb.hit")
    primary = resp.get("bdb.primary")
    logger.info(f"BindingDB reported hits: {hits}, primary UniProt: {primary}")

    rows = resp.get("bdb.affinities", [])
    if not isinstance(rows, list):
        raise TypeError("bdb.affinities is not a list")

    df = pd.DataFrame(rows)
    logger.info(f"Rows in bdb.affinities: {len(df)}")
    return df

# clean smiles
def clean_smiles(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip()
    # Remove CXSMILES / coordinate blocks like: "SMILES |...|"
    if "|" in s:
        s = s.split("|", 1)[0].strip()
    # Sometimes there is extra whitespace before the pipe
    return s

def run_bindingdb_pipeline(cfg: BDBConfig) -> pd.DataFrame:
    logger = setup_logger(cfg.log_path)
    logger.info("=== BindingDB pipeline start ===")
    logger.info(f"UniProt: {cfg.uniprot}")
    logger.info(f"Assumed IC50 unit: {cfg.ic50_unit}")

    df = fetch_bindingdb_affinities(cfg.uniprot, logger, cfg)

    # Expected columns from your snippet:
    # bdb.monomerid, bdb.smile, bdb.affinity_type, bdb.affinity
    required = ["bdb.monomerid", "bdb.smile", "bdb.affinity_type", "bdb.affinity"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"Missing expected columns: {missing}. Got: {list(df.columns)}")

    # Clean whitespace / types
    logger.info("[Step] Clean fields")
    before = len(df)
    df = df[df["bdb.affinity_type"]=="IC50"]
    logger.info(f"Number of IC50 data: {len(df)}")
    df["SMILES"] = df["bdb.smile"].apply(clean_smiles)

    # Clean uncertain data
    df = (
    df[~df["bdb.affinity"].astype(str).str.contains("<|>", regex=True)]
      .assign(IC50=lambda d: d["bdb.affinity"].str.extract(r"(\d+\.?\d*)").astype(float))
      .dropna(subset=["IC50"]))
    
    # Numeric affinity
    logger.info("[Step] Parse numeric IC50 and compute pIC50")
    df["pIC50"] = compute_pic50(df["IC50"], cfg.ic50_unit)
    df["pIC50"] = df["pIC50"].round(2)
    df = df.drop(columns=["bdb.smile", "bdb.affinity_type", "bdb.affinity", "IC50"])
    logger.info(f"Rows after basic cleanup: {before} -> {len(df)}")

    logger.info(
        f"pIC50 stats: n={len(df)}, min={df['pIC50'].min():.2f}, "
        f"median={df['pIC50'].median():.2f}, max={df['pIC50'].max():.2f}"
    )

    # Deduplicate (median pIC50 per SMILES)
    before = len(df)
    print(df)
    out = (df
           .groupby("SMILES", as_index=False)["pIC50"]
           .median())
    print(out)
    logger.info(f"Unique Monomer ID: {before} rows -> {len(out)} ID")

    if cfg.out_csv:
        cfg.out_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(cfg.out_csv, index=False)
        logger.info(f"Wrote: {cfg.out_csv}")

    logger.info("=== BindingDB pipeline done ===")
    return out


# -----------------------------
# CLI
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="BindingDB: UniProt -> IC50 -> SMILES + median pIC50")
    parser.add_argument("--uniprot", required=True, help="UniProt accession (e.g., P27815)")
    parser.add_argument("--out_csv", default=None, help="Output CSV path")
    parser.add_argument("--log_path", default="bindingdb_pipeline.log", help="Log file path")
    parser.add_argument("--ic50_unit", default="nM", choices=["uM", "nM", "M"], help="Unit of bdb.affinity values")
    parser.add_argument("--pause_sec", type=float, default=0.25, help="Delay between requests")

    args = parser.parse_args()

    cfg = BDBConfig(
        uniprot=args.uniprot,
        out_csv=Path(args.out_csv) if args.out_csv else None,
        log_path=args.log_path,
        ic50_unit=args.ic50_unit,
        pause_sec=args.pause_sec,
    )
    run_bindingdb_pipeline(cfg)


if __name__ == "__main__":
    main()
