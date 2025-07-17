# Functional connectomics reveals general wiring rule in mouse visual cortex

[![DOI](https://zenodo.org/badge/DOI/10.1038/s41586-025-08840-3.svg)](https://doi.org/10.1038/s41586-025-08840-3)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/cajal/microns-funconn-2025)](https://github.com/cajal/microns-funconn-2025/releases)
[![GitHub license](https://img.shields.io/github/license/cajal/microns-funconn-2025)](https://github.com/cajal/microns-funconn-2025/blob/main/LICENSE)

This repository contains the code for the paper "Functional connectomics reveals general wiring rule in mouse visual cortex".

[https://doi.org/10.1038/s41586-025-08840-3](https://doi.org/10.1038/s41586-025-08840-3).

## Installation

You can set up this project in several ways:

### Option 1: Using Docker (Recommended)

The easiest way to get started is to use our pre-built Docker container which includes all required dependencies.

#### Prerequisites
- [Docker](https://docs.docker.com/get-docker/)
- [Docker Compose](https://docs.docker.com/compose/install/) (included in Docker Desktop for Mac/Windows)

#### Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/cajal/microns-funconn-2025.git
   cd microns-funconn-2025
   ```

2. Start the Docker container:
   ```bash
   ./run-docker-container.sh
   ```

3. Access Jupyter Lab:
   Open your browser and navigate to `http://localhost:8888`
   Note: It may take a few seconds for the content to become accessible.


### Option 2: Install with pip

#### Prerequisites
- Python 3.8 or higher
- R 4.0.0 or higher (for statistical analysis)

#### Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/cajal/microns-funconn-2025.git
   cd microns-funconn-2025
   ```

2. Install the package and its Python dependencies:
   ```bash
   # Install in development mode
   pip install -e .
   ```

3. Install required R packages:
   ```bash
   # Run the R package setup script
   Rscript setup_r_packages.R
   ```

   This script will install the following R packages:
   - glmmTMB
   - tidyverse
   - broom.mixed
   - emmeans
   - performance
   - DHARMa


## Reproduce the figures

The intermediate results files are already included in the `results` folder. Notebooks to load these files and reproduce the figures in the paper are under the `figures` folder.

1. Navigate to the `figures` directory
2. Run the Jupyter notebooks:
   - `like2like.ipynb` - Figures related to the like-to-like connectivity analysis
   - `common_input.ipynb` - Figures related to the common input analysis

To reproduce the intermediate results, you can run the following scripts:
- `funconnect/compute/like2like.py` - Script to generate the like-to-like connectivity analysis results
- `funconnect/compute/common_inputs.py` - Script to generate the common input analysis results

To run the scripts, open a terminal inside the `./funconnect/compute/` directory and run:
```bash
python3 ./common_inputs.py
python3 ./like2like.py
```

Intermediate results are stored in `funconnect/results` and should match the results in the `results` folder.

## Reproduce methods

Notebooks to demonstrate some of the methods used in the paper are in the `methods` folder. Currently available methods:
- proximities - `compute_proximities.ipynb`

## Data Availability

To access the datasets analyzed in this study, please see the [data availability](https://www.nature.com/articles/s41586-025-08840-3#data-availability) section of the manuscript.

They are also downloaded inside the Docker container at `/data`.

## Citation

If you find this repository useful, please cite using this BibTeX:

```bibtex
@article{ding2025functional,
  title={Functional connectomics reveals general wiring rule in mouse visual cortex},
  author={Ding, Zhuokun and Fahey, Paul G and Papadopoulos, Stelios and Wang, Eric Y and Celii, Brendan and Papadopoulos, Christos and Chang, Andersen and Kunin, Alexander B and Tran, Dat and Fu, Jiakun and others},
  journal={Nature},
  volume={640},
  number={8058},
  pages={459--469},
  year={2025},
  publisher={Nature Publishing Group UK London}
}
```