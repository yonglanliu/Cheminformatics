# Cheminformatics Data Pipeline
This repository provides a streamlined workflow for fetching, cleaning, and processing bioactivity and molecular structure data directly from the ChEMBL database.

##🚀 Getting Started

1. Environment Setup
Ensure you have Conda installed. Create the environment from the provided environment.yml file:

## Conda Environment setup and activate environment

``bash
# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate chem
``


###🧬 Data Acquisition Pipeline

The pipeline automates the retrieval of bioactivity data, filters for IC50 values, handles missing data, and converts values to $pIC_{50}$ ($-\log_{10}(IC_{50})$).

Configuration
Before running the pipeline, update the variables in <code>./bash/fetch_chembl_data.sh</code> to match your target of interest:


``bash
TARGET_CHEMBL_ID=" "  # The ChEMBL ID of your protein target (e.g., CHEMBL203).

BIOACTIVITY_OUT_PATH="bioactivity_chembl.csv"  # (Optional) Filename for the raw bioactivity data download.
FINAL_OUT_PATH="smiles_plus_pIC50.csv"   # Path for the cleaned CSV containing SMILES and $pIC_{50}$.
LOG_PATH="chembl_pipeline.log"  # Location for the process logs.

RUN_SCRIPT="./scripts/fetch_chembl_data.py"
``

``bash
source ./bash/fetch_chembl_data.sh
``

### 🛠 Features

Data Cleaning: Automatic removal of duplicates and NaN values.Standardization: Filters specifically for IC50 bioactivity types.Transformation: Automated conversion of molar IC50 concentrations to $pIC_{50}$ for better statistical distribution.Logging: Detailed step-by-step tracking of the retrieval and cleaning process.
