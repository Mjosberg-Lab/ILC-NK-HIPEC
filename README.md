[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18618925.svg)](https://doi.org/10.5281/zenodo.18618925)
# ILC-NK-HIPEC Manuscript Code

This repository contains the analysis code used to generate the figures for our manuscript "Tumor-infiltrating immature innate lymphoid cells in colorectal cancer are biased towards tissue-resident NK cell/ILC1 differentiation" (Marchalot et al.)

## Data Availability

The analyzed Seurat object `ILC-NK` required to run most of the code is available on GEO with accession number **GSE302045**.

## Repository Contents

This repository contains R and Python scripts organized by the figures they generate:

### R Scripts
- `00_setup.R` - Package installation and loading
- `Figure 1.R` - Plotting code for Figure 1
- `Figure 2.R` - Plotting code for Figure 2
- `Figure 3A-D.R` - Plotting code for Figure 3 panels A-D
- `Figure 4.R` - Plotting code for Figure 4
- `Figure 5A-B.R` - Plotting code for Figure 5 panels A-B
- `Figure 5C-D.R` - Plotting code for Figure 5 panels C-D

### Python Notebooks
- `Figure 2E-F.ipynb` - RNA velocity analysis for Figure 2 panels E-F
- `Figure 3D.ipynb` - RNA velocity analysis for Figure 3 panel D
- `Figure 6H-J.ipynb` - Analysis code for Figure 6 panels H-J

## Requirements

### R Dependencies
See `00_setup.R` for the complete list of required R packages including:
- Seurat (single-cell analysis)
- harmony (batch correction)
- monocle3 (trajectory analysis)
- And many others listed in the setup file

### Python Dependencies
For the Jupyter notebooks (velocity analysis):
- scanpy
- scvelo
- numpy
- pandas
- matplotlib

## Usage

1. Download the analyzed data object from GEO (GSE302045)
2. Install required packages using `00_setup.R`
3. Run the individual figure scripts with the appropriate data paths

## Note on RNA Velocity Analysis

The RNA velocity analyses (Figure 2E-F and Figure 3D) require loom files that contain spliced/unspliced RNA counts. Due to patient privacy concerns, these raw loom files cannot be shared as they contain sensitive information. The code is provided for transparency and reproducibility purposes, showing the exact methods used in the manuscript.

## Citation

If you use this code, please cite our manuscript: [Citation to be added upon publication]

## Contact

For questions about the code or data, please open an issue in this repository.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
