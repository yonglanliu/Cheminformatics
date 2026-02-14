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
PROTEIN_NAME=""

PARENT_PATH=""

CHEMBL_DATA_PATH="$PARENT_PATH/data/${PROTEIN_NAME}/smiles_plus_pIC50_chembl.csv"
PUBCHEM_DATA_PATH="$PARENT_PATH/data/${PROTEIN_NAME}/smiles_plus_pIC50_pubchem.csv"
BINDINGDB_DATA_PATH="$PARENT_PATH/data/${PROTEIN_NAME}/smiles_plus_pIC50_BindingDB.csv"
OUT_PATH="$PARENT_PATH/data/${PROTEIN_NAME}/out.csv"
OUT_FINAL="$PARENT_PATH/data/${PROTEIN_NAME}/combine_database.csv"
LOG_PATH="$PARENT_PATH/data/${PROTEIN_NAME}/combine_database.log"
OUT_OVERLAP="$PARENT_PATH/data/${PROTEIN_NAME}/overlap_by_source.csv"

RUN_SCRIPT="./scripts/combine_database.py"

# ----------------------------
# Run pipeline
# ----------------------------
python "$RUN_SCRIPT" \
    --chembl "$CHEMBL_DATA_PATH" \
    --pubchem "$PUBCHEM_DATA_PATH" \
    --bindingdb "$BINDINGDB_DATA_PATH" \
    --out_combined "$OUT_PATH" \
    --out_final "$OUT_FINAL" \
    --out_overlap "$OUT_OVERLAP"\
    --log_path "$LOG_PATH"

echo "Pipeline finished successfully."
