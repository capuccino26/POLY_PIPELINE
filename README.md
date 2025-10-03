# POLY_PIPELINE

> A data analysis pipeline for STOmics data tailored to **polyploid** organisms.

---

## 📘 Overview

**POLY_PIPELINE** is a comprehensive bioinformatics pipeline designed to process and analyze **STOmics** data specifically tailored for **polyploid** organisms. It automates and simplifies several crucial steps of spatial transcriptomics data analysis, ensuring reproducibility, flexibility, and scalability for handling large, complex polyploid genomes. This pipeline leverages a modular structure and relies heavily on the **Stereopy**.

---

## 🤝 Contributing

Contributions are welcome! Please fork the repository and submit a pull request.
See the [`CONTRIBUTING.md`](CONTRIBUTING.md) for details.

---

## 🛠️ Usage

### Initial Setup

The pipeline is designed to be executed from the **main folder** (`POLY_PIPELINE/`) and requires a dedicated Conda environment.

1.  **Set Environment:**
    * Run `bin/0_SET_ENV.sh` to activate the **`st`** Conda environment, which is essential due to the strict package dependencies of **Stereopy**.
2.  **Test Environment:**
    * Run `bin/1_TEST_ENV.sh` to verify the installation.

---

### Main Pipeline (Cluster Execution - SGE)

The core scripts are optimized for **SGE cluster** execution. They use **relative paths** and **must be run from the `POLY_PIPELINE` main directory.**

####  Data Input
The input file must be the `.gef` file (post-processed by the SAW pipeline). It **must be placed** in the `INPUT/datasets/` folder.
> **IMPORTANT:** Place **only one** `.gef` file in the `datasets` folder.
It is possible **(BUT NOT REQUIRED)** to generate differential analysis for a list of genes of interest, generate the file `INPUT/datasets/interest_genes.txt` following the structure:
```markdown
gene_name,Gene_ID_1,Gene_ID_2,Gene_ID_3,Gene_ID_4,Gene_ID_5
FLORAL_MERISTEM,AT5G08570,LOC107775591,Nicotiana_T001,LOC107775592
STRESS_HEAT,AT1G53540,LOC107817066,GmHIS4_A01,LOC107817067,AT2G41090
AUXIN_RESPONSE,AT3G15540,LOC107833544,Os02g0602300,AT4G20560
CELL_CYCLE,LOC107769919,AT1G44110,LOC107769920,AT3G53210
APICAL_DOMINANCE,LOC107802111,AT2G44320,LOC107802112
DEFENSE_MECH,AT5G41220,LOC107764120,Solyc01g099710
GIBBERELLIN_SYN,LOC107823450,AT1G05030,LOC107823451,AT3G44360
```
Each line represents one gene of interest starting with the identification of the gene followed by all correponding IDs (The IDs must match the mapping reference used in the generation of the .gef file).

| Step | Script | Description |
| :--- | :--- | :--- |
| [**1. Analysis**](bin/2_DOC_ANALYSIS.sh) | `bin/2_DOC_ANALYSIS.sh` | Complete analysis following the Stereopy documentation (generates the `stereopy_ultimate_analysis.py` script). |
| [**2. Annotation**](bin/3_ANNOTATION.sh) | `bin/3_ANNOTATION.sh` | Annotates results from Step 1. |
| [**3. Validation**](bin/4_VALIDATEH5.sh) | `bin/4_VALIDATEH5.sh` | Validates the output `.h5ad` file against the original file (generates `data_validation_analysis.py`). **Note:** "Failed analysis" due to gene filtering is **normal**; this step serves as a reference. |
| [**4. Packaging**](bin/5_ZIP_RESULTS.sh) | `bin/5_ZIP_RESULTS.sh` | Zips the result folders for local transfer and downstream analysis. |

#### Cluster Execution Example (SGE)

The scripts are submitted with explicit Miniconda and parameter variables (`qsub -v`).

* [**Analysis Script:**](bin/2_DOC_ANALYSIS.sh)
    ```bash
    qsub -v ST_PYTHON="/home/user/.conda/envs/st/bin/python",MIN_COUNTS=50,MIN_GENES=5,PCT_COUNTS_MT=100,N_PCS=30 bin/2_DOC_ANALYSIS.sh
    ```
    | Variable | Description |
    | :--- | :--- |
    | `ST_PYTHON` | Path to the python executable inside the st environment. |
    | `MIN_COUNTS` | Minimum number of counts per cell. |
    | `MIN_GENES` | Minimum number of genes per cell. |
    | `PCT_COUNTS_MT` | Acceptable percentage of mitochondrial genes. |
    | `N_PCS` | Number of principal components. This step can be inproved after first run. Ceck the Elbow Plot (RESULTS/results_ultimate/plots/qc/pca_elbow_enhanced.png) and insert the value of the elbow as N_PCS |

#### Local Execution Example
> **IMPORTANT:** This analysis requires high computational resources and are not recommended to be run locally.
To run the main analysis locally using the `bash` wrapper and your specific Conda path:

```bash
ST_PYTHON='/home/user/.conda/envs/st/bin/python' MIN_COUNTS=50 MIN_GENES=5 PCT_COUNTS_MT=100 N_PCS=30 bash bin/2_DOC_ANALYSIS.sh
```

---

## 📂 Project Structure
```markdown

POLY_PIPELINE/
├── CONTRIBUTING.md
├── INPUT/
│   └── datasets/
│       └── <place .gef file here>
├── LICENSE
├── RESULTS/
│   └── <results folders and .zip files generated>
├── README.md
└── bin/
    ├── 0_SET_ENV.sh
    ├── 1_TEST_ENV.sh
    ├── 2_DOC_ANALYSIS.sh
    ├── 3_ANNOTATION.sh
    ├── 4_VALIDATEH5.sh
    └── 5_ZIP_RESULTS.sh
```

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).

---
