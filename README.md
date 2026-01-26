# POLY_PIPELINE

> A data analysis pipeline for STOmics data tailored to **polyploid** organisms.

---

## Overview

**POLY_PIPELINE** is a comprehensive bioinformatics pipeline designed to process and analyze **STOmics** data specifically tailored for **polyploid** organisms. It automates and simplifies several crucial steps of spatial transcriptomics data analysis, ensuring reproducibility, flexibility, and scalability for handling large, complex polyploid genomes. This pipeline leverages a modular structure and relies heavily on the [**Stereopy package**](https://github.com/STOmics/Stereopy).

---

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request.
See the [`CONTRIBUTING.md`](CONTRIBUTING.md) for details.

---

## Usage

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
The input file must be the [`.gef`](https://www.processon.com/view/link/610cc49c7d9c087bbd1ab7ab#map) file (post-processed by the [**SAW pipeline**](https://github.com/STOmics/SAW)). It **must be placed** in the `INPUT/datasets/` folder.
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
| [**Analysis AND Annotation**](bin/2_COMPLETE_ANALYSIS.sh) | `bin/2_COMPLETE_ANALYSIS.sh` | Complete analysis following the Stereopy documentation (generates the `stereopy_ultimate_analysis.py` script). |

#### Cluster Execution Example (SGE)

The scripts are submitted with explicit Miniconda and parameter variables (`qsub -v`).

* [**Analysis Script:**](bin/2_COMPLETE_ANALYSIS.sh)
    ```bash
    qsub -v ST_PYTHON="/home/user/.conda/envs/st/bin/python",MIN_COUNTS=50,MIN_GENES=5,PCT_COUNTS_MT=100,N_PCS=30 bin/2_COMPLETE_ANALYSIS.sh
    ```
    **Optional: Add coordinate filtering to analyze specific regions**
    ```bash
    qsub -v ST_PYTHON="/home/user/.conda/envs/st/bin/python",MIN_COUNTS=50,MIN_GENES=5,PCT_COUNTS_MT=100,N_PCS=30,MIN_X=7176,MAX_X=16425,MIN_Y=5300,MAX_Y=12200 bin/2_COMPLETE_ANALYSIS.sh
    ```
    | Variable | Description |
    | :--- | :--- |
    | `ST_PYTHON` | Path to the python executable inside the st environment. |
    | `MIN_COUNTS` | Minimum number of counts per cell. |
    | `MIN_GENES` | Minimum number of genes per cell. |
    | `PCT_COUNTS_MT` | Acceptable percentage of mitochondrial genes. |
    | `N_PCS` | Number of principal components. This step can be inproved after first run. Check the Elbow Plot (RESULTS/results_ultimate/plots/qc/pca_elbow_enhanced.png) and insert the value of the elbow as N_PCS |
    | `MIN_X` | *(Optional)* Minimum X coordinate for spatial filtering. |
    | `MAX_X` | *(Optional)* Maximum X coordinate for spatial filtering. |
    | `MIN_Y` | *(Optional)* Minimum Y coordinate for spatial filtering. |
    | `MAX_Y` | *(Optional)* Maximum Y coordinate for spatial filtering. |

* The variables are not required, the script can run with defaults and the entire tissue area.
* If coordinate filtering is required (MIN_X, MAX_X, MIN_Y, MAX_Y), all coordinate parameters must be provided together.
#### Local Execution Example
> **IMPORTANT:** This analysis requires high computational resources and are not recommended to be run locally.
To run the main analysis locally using the `bash` wrapper and your specific Conda path:

```bash
ST_PYTHON='/home/user/.conda/envs/st/bin/python' MIN_COUNTS=50 MIN_GENES=5 PCT_COUNTS_MT=100 N_PCS=30 MIN_X=7176 MAX_X=16425 MIN_Y=5300 MAX_Y=12200 bash bin/2_COMPLETE_ANALYSIS.sh
```

---

### Miscellaneous Pipeline (Secondary analysis)

| Step | Script | Description | Usage | Observations |
| :--- | :--- | :--- | :--- | :--- |
| [**Plotting interest genes over sample**](bin/MISC_01_PLOT_CELL.py) | `bin/MISC_01_PLOT_CELL.py` | Script for generating plots with gene of interest expression overlay over sample. | python bin/MISC_01_PLOT_CELL.py -i INT_GENES.txt -o CLUSTER10 -t 10 --min_x 7176 --max_x 16425 --min_y 5300 --max_y 12200 | Only variable -i (--interest) is required, the others are optional; Check example file below |

* For script [**Plotting interest genes over sample**](bin/MISC_01_PLOT_CELL.py) the gene file list should be as below, with one gene/loc per line:
```markdown
LOC123047130
LOC123091185
LOC123147796
LOC123112488
LOC123126477
LOC123045745
LOC123129052
```
* If coordinate filtering is required (min_y, max_y, min_x, max_x), all coordinate parameters must be provided together.
---

## Project Structure
```markdown

POLY_PIPELINE/
├── CONTRIBUTING.md
├── INPUT/
│   └── datasets/
│       ├── <place .gef file here>
│       └── <(OPTIONAL) place interest_genes.txt here>
├── LICENSE
├── RESULTS/
│   └── <results folders and .zip files generated>
├── README.md
└── bin/
    ├── 0_SET_ENV.sh
    ├── 1_TEST_ENV.sh
    ├── 2_COMPLETE_ANALYSIS.sh
    └── MISC_01_PLOT_CELL.py

```

---

## License

This project is licensed under the [MIT License](LICENSE).

---
