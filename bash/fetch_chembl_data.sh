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
TARGET_CHEMBL_ID=""

BIOACTIVITY_OUT_PATH=""
FINAL_OUT_PATH=""
LOG_PATH=""

RUN_SCRIPT=""

# ----------------------------
# Run pipeline
# ----------------------------
python "$RUN_SCRIPT" \
    --target_chembl_id "$TARGET_CHEMBL_ID" \
    --out_final_csv "$FINAL_OUT_PATH" \
    --out_bio_csv "$BIOACTIVITY_OUT_PATH" \
    --log_path "$LOG_PATH"

echo "Pipeline finished successfully."
