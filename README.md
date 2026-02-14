# Cheminformatics Data Pipeline
This repository provides a streamlined workflow for fetching, cleaning, and processing bioactivity and molecular structure data directly from different database.

## Getting Started

1. Environment Setup
Ensure you have Conda installed. Create the environment from the provided environment.yml file:

## Conda Environment setup and activate environment

```bash
# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate chem
```


### Data Acquisition Pipeline

The pipeline automates the retrieval of bioactivity data, filters for IC50 values, handles missing data, and converts values to <code>pIC50 = -log10(IC50_Molar)</code>.

### Configuration
Before running the pipeline, update the variables in <code>./bash/fetch_chembl_data.sh</code> or <code>./bash/fetch_pubchem_data.sh</code> to match your target of interest:

```bash
# --- Configuration Settings ---

TARGET_CHEMBL_ID=""  # For ChEMBL data --- Target ChEMBL ID (e.g., CHEMBL203 for EGFR)
TARGET_UNIPROT_ID="" # For PubChem and BindingDB data

# Output paths
BIOACTIVITY_OUT_PATH="bioactivity.csv"        # Raw data storage
FINAL_OUT_PATH="smiles_plus_pIC50.csv"        # Processed molecule_chembl_id/CID + SMILES + pIC50. Note: output of Binding DB data just have SMILES + pIC50
LOG_PATH="pipeline.log"                       # Execution logs

# Script location
RUN_SCRIPT="./scripts/fetch_chembl_data.py"  # For ChEMBL dataset
RUN_SCRIPT="./scripts/fetch_pubchem_data.py"  # For PubChem dataset
RUN_SCRIPT="./scripts/fetch_BindingDB_data.py"  # For BindingDB dataset
```

Then run:


```bash
source ./bash/fetch_chembl_data.sh # To fetch chembl data
source ./bash/fetch_pubchem_data.sh # To fetch pubchem data
source ./bash/fetch_BindingDB_data.sh # To fetch BindingDB data
```

### Features

* **Data Cleaning**: Automatic removal of duplicates and NaN values. For the multiple data with the same molecule_chembl_id, we use the median value of pIC50
* **Standardization**: Filters specifically for IC50 bioactivity types.
* **Transformation**: Automated conversion of molar IC50 concentrations to $pIC_{50}$ for better statistical distribution.
* **Logging**: Detailed step-by-step tracking of the retrieval and cleaning process.


## Combine Databases
**Workflow:**

1. Clean SMILES (remove |...|)
2. Convert SMILES → RDKit Mol
3. Canonicalize (optional but good)
4. Generate InChIKey
5. Merge on InChIKey


**Note**: For a molecule, InChIKey is unique, but SMILES may change.

✅ InChIKey is intended to uniquely represent a chemical structure (after standardization).

#### ⚠️ SMILES can vary for the same molecule.

* SMILES is a string encoding of a graph, and:
* There are many valid SMILES for the same molecule.
* Different databases use different canonicalization rules.
* Atom ordering can change.
* Explicit vs implicit hydrogens can change.
* Aromatic vs Kekulé form can change.

Example (same molecule, caffeine):

```bash
CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Cn1cnc2n(C)c(=O)n(C)c(=O)c12
```

Different SMILES, same molecule. So merging datasets on SMILES is risky unless you canonicalize first.

#### InChIKey is:

* Derived from the standardized InChI
* Canonical
* Structure-based
* Database-independent
* Designed for cross-database matching

For the same standardized structure:

```bash
ZFXYFBGIUFBOJW-UHFFFAOYSA-N
```
