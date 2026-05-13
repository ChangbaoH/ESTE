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

#### Step 2: Estimate Error Rates (epsilon)
```r
# Fast epsilon estimation
epsilon_result <- estimate_Epsilon(
  pat = as.integer(genotype_data),
  isF = isF,
  isCE = isCE,
  eps = eps,
  setD = setD,
  eventD = eventD,
  threshold = 0.05,
  thrds = 4  # Number of parallel threads
)
```

#### Step 3: Infer Partial Order (Poset)
```r
# Find poset with fine-tuning
poset_result <- find_Poset(
  pat = as.integer(genotype_data),
  isF = isF,
  isCE = isCE,
  eps = epsilon_result,
  setD = setD,
  eventD = eventD,
  Fine_Tune_Num = 2L,
  threshold = 0.0001,
  thrds = 4
)
```

#### Step 4: Estimate Rate Parameters (lambda)
```r
# Estimate lambda using EM algorithm
lambda_result <- estimate_Lambda(
  pat = as.integer(genotype_data),
  poset = poset_result,
  isF = isF,
  eps = epsilon_result,
  setD = setD,
  eventD = eventD,
  sampling = "add-remove",
  maxIter = 100L,
  thrds = 4
)
```

#### Step 5: Sample Event Timing (Optional)
```r
# MH sampling for event timing estimation
timing_result <- sample_Age_T(
  data = genotype_data,
  alpha = 1.0,
  beta = 1.0,
  lambda = lambda_result$lambda,
  maxSampleIter = 1000,
  rateE1 = 0.01,
  thrds = 4,
  seed = 12345
)
```

### Advanced Usage

#### Multiple Datasets and Event Sets
```r
# For multi-dataset analysis
epsilon_multi <- estimate_Epsilon_ForMulti(
  pat = as.integer(genotype_data),
  isF = isF_matrix,
  isCE = isCE_matrix,
  eps = eps_matrix,
  setD = setD_matrix,
  eventD = eventD_matrix,
  multi_thrds = 4,
  thrds = 2
)

# Voting for consensus poset
poset_vote <- find_Poset_ForVote(
  pat = as.integer(genotype_data),
  isF = isF,
  eps = epsilon_multi,
  setD = setD,
  eventD = eventD,
  vote_size = 5L,
  vote_threshold = 0.5,
  thrds = 4
)
```

#### Compatibility Check
```r
# Check if genotypes are compatible with poset
compatibility <- is_Compatible(
  genotype = genotype_matrix,
  poset = poset_result
)
```

---

## 5. Reproduction Instructions

To reproduce the results from the associated research paper:

### Step 1: Download Full Dataset
```bash
# Download the complete ICGC dataset (if not already present)
# The dataset should be placed in data/ICGC/ directory
```

### Step 2: Run Analysis Scripts
```r
# Load required packages
library(este)
library(parallel)
library(data.table)
library(dplyr)

# Run pan-cancer CBN analysis
source("example/panCancerCBNandMH-Sampling.R")

# Run fitness inference
source("example/panCancerImmuneFitnessInference.R")

# Run single-cell RNA analysis
source("example/scRNA_Analysis.R")
```

### Step 3: Parameter Settings

The following parameters were used in the original analysis:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Fine_Tune_Num` | 2 | Number of fine-tuning iterations |
| `threshold` | 0.0001 | Poset inference threshold |
| `threshold2` | 0.01 | Tolerance threshold |
| `sampling` | "add-remove" | Hidden genotype sampling method |
| `maxIter` | 100 | EM algorithm iterations |
| `L` | 100 | Number of E-step samples |
| `tol` | 0.001 | Convergence tolerance |

### Step 4: Output Interpretation

**Partial Order Matrix**:
- `poset[i,j] = 1` indicates event i must occur before event j
- Diagonal elements are always 0
- The matrix represents a directed acyclic graph (DAG)

**Lambda Values**:
- Higher lambda = faster event occurrence
- Lambda values are relative rates

**Timing Estimates**:
- Estimated time of each event occurrence
- Relative timing between events indicates progression order

---

## 6. Package Structure

```
este/
├── R/                    # R source code
│   ├── este.R            # Main functions (epsilon, poset, lambda)
│   ├── baseFunction.R    # Utility functions (topological sort, etc.)
│   ├── baselineTest.R    # Simulation and testing functions
│   └── RcppExports.R     # Rcpp function exports
├── src/                  # C++ source code
│   ├── este.cpp          # Main C++ implementation
│   ├── ct-cbn.h          # CT-CBN algorithm header
│   ├── gmm.h             # Gaussian Mixture Model header
│   ├── lambdaCal.h       # Lambda calculation header
│   ├── interfaceFun.h    # R-C++ interface functions
│   └── rng_utils.hpp     # Random number generation utilities
├── data/                 # Example datasets
│   └── cancer_data/      # Cancer genotype matrices
├── example/              # Example scripts
├── man/                  # R documentation files
├── DESCRIPTION           # Package metadata
├── NAMESPACE             # Export declarations
└── README.md             # This file
```

---

## 7. API Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `estimate_Epsilon()` | Estimate error rates for hidden CBN model |
| `estimate_Epsilon_based_on_Poset_and_Lambda()` | Epsilon estimation with known poset/lambda |
| `estimate_Epsilon_ForMulti()` | Multi-dataset epsilon estimation |
| `find_Poset()` | Infer partial order from genotype data |
| `find_Poset_ForVote()` | Consensus poset via voting |
| `estimate_Lambda()` | Estimate event rate parameters |
| `estimate_Lambda_ForMulti()` | Multi-dataset lambda estimation |
| `is_Compatible()` | Check genotype-poset compatibility |
| `sample_Age_T()` | MH sampling for event timing |
| `sample_Age_TLW()` | Lightweight timing sampling |

### Utility Functions

| Function | Description |
|----------|-------------|
| `topological_Sort()` | Perform topological sort on poset |
| `random_Poset()` | Generate random poset structure |
| `rateTimeHelp()` | Convert lambda to time bounds |
| `find_most_Compatible_Genotype_by_Flipping()` | Find compatible genotypes |
| `donor_Pair_Genotype_Filter()` | Filter donor pair genotypes |
| `GammaCluster()` | Gaussian Mixture Model clustering |

### Data Generation Functions

| Function | Description |
|----------|-------------|
| `simulation_Data_Generate()` | Generate simulated genotype data |
| `simulation_Time_Data_Generate()` | Generate simulated timing data |

---

## 8. Citation

If you use ESTE in your research, please cite:

```
Hu, C. (2024). ESTE: Estimate the Sequence and Timing of Events. 
R package version 0.2.0. https://github.com/yourusername/este
```

---

   
## Support
For questions and support, contact: Hu Changbao <1437894182@qq.com>

## License
GPL (≥ 2)
