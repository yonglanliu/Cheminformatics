#!/usr/bin/env bash

# ----------------------------
# Load conda
# ----------------------------
source ~/bin/myconda
conda activate chem

echo "Conda environment activated: $(which python)"

# ----------------------------
# Variables
# ----------------------------
TARGET_UNIPROT_ID="" # Use uniprot ID

FINAL_OUT_PATH="smiles_plus_pIC50_BindingDB.csv"
LOG_PATH="BindingDB_pipeline.log"

RUN_SCRIPT="./scripts/fetch_BindingDB_data.py"

# ----------------------------
# Run pipeline
# ----------------------------
python "$RUN_SCRIPT" \
    --uniprot "$TARGET_UNIPROT_ID" \
    --out_csv "$FINAL_OUT_PATH" \
    --log_path "$LOG_PATH"

echo "Pipeline finished successfully."
