# ESTE
ESTE is an R package that implements the improved CBN algorithm to estimate the sequence of events and MH-Sampling to estimate the timing of events. This package is specifically designed for cancer genomic data analysis, providing tools for inferring the temporal order of genetic alterations in cancer progression.

## 1. System Requirements
### Software Dependencies
- **R**: Version 4.4.1
- **R Packages**:
  - Rcpp (1.0.12)
  - relations (0.6-13)
  - BH (1.84.0-0)
  - RcppEigen (0.3.4.0.0)
### Operating Systems
- Linux (Ubuntu 24.04 LTS)
### Tested Versions
- R 4.4.1 on Ubuntu 24.04
### Hardware Requirements
- No special hardware requirements
  
## 2. Installation Guide
### Installation Instructions
1. **Install R and required dependencies** :
2. **Install ESTE package from source**:
3. **Load the package**:

## 3. Example
### Instructions to Run Demo
The package includes example data and scripts in the `example/` and `data/` directory.
### Dataset
The real dataset is provided in `data/cancer_data/` containing genotype matrices for various cancer types:
- Breast-AdenoCA
- CNS-GBM  
- ColoRect-AdenoCA
- Liver-HCC
- Lung-AdenoCA
- Lung-SCC
- Prost-AdenoCA
- Skin-Melanoma
- Uterus-AdenoCA
  
## 4. Instructions for Use
### Data Preparation
ESTE requires genotype data in a specific format:
```r
# Genotype matrix format:
# Rows: samples/patients
# Columns: genetic events (binary: 1 = observed, 0 = not observed)
# Last column: dataset identifier (optional)

# Example:
genotype_data <- read.csv("your_genotype_data.csv", row.names = 1)
```

### Basic Usage Workflow

#### Step 1: Prepare Input Matrices
```r
# Define dataset and event set configurations
setD <- matrix(c(0, nrow(genotype_data)-1), nrow = 1)  # Dataset range
eventD <- matrix(c(0, ncol(genotype_data)-2), nrow = 1)  # Event range
isF <- matrix(1L, nrow = 1, ncol = 1)  # Data filled indicator
isCE <- matrix(1L, nrow = 1, ncol = 1)  # Calculate epsilon indicator
eps <- matrix(0.05, nrow = 1, ncol = 1)  # Initial epsilon estimate
```

### Running on Your Data
1. **Prepare your genotype data**:
   - Format as CSV with events as columns and samples as rows
   - Binary encoding: 1 = event observed, 0 = not observed
   - Include dataset and sample identifiers
2. **Basic usage**:
   please see example/baselineTestExample_este.R
   
## Support
For questions and support, contact: Hu Changbao <1437894182@qq.com>

## License
GPL (≥ 2)
