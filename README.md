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

#### Please check example\baselineTestExample_este.R

### Generate Simulated Data

ESTE provides functions to generate simulated genotype data for testing and validation purposes.

```r
# Load the package
library(este)

# Generate simulated data
# Parameters:
#   numEventArray: Array describing the number of events in each event type
#   numSampleArray: Array describing the number of samples in each dataset
#   baseEpsilon: Base error rate
#   lambdaSampling: Sampling rate parameter
#   lambdaSamplingScaling: Scaling factor for lambda sampling
#   epsilonSamplingRate: Rate for epsilon variation
#   graphDensity: Density of the poset graph

# Example: Generate data with 2 event types (8 events each) and 2 datasets (800 samples each)
numEventKind <- c(2)  # 2 event types
numEventPerKind <- c(8)  # 8 events per event types
numSet <- c(2)  # 2 datasets
numSamplePerSet <- 800  # 800 samples per dataset

sim_data <- simulation_Data_Generate(
  numEventArray = rep(numEventPerKind[1], numEventKind),  # Event array
  numSampleArray = rep(numSamplePerSet, numSet),                     # Sample array
  baseEpsilon = 0.05,           # Base error rate
  lambdaSampling = 1,           # Sampling rate
  lambdaSamplingScaling = 3,    # Lambda scaling factor
  epsilonSamplingRate = 0.5,    # Epsilon variation rate
  graphDensity = 0.2            # Poset graph density
)
```

**Output Structure**:

- `obs_events`: Matrix of observed genotypes with noise (rows = samples, cols = events)

- `hidden_genotypes`: True underlying genotypes without noise

- `poset`: Ground truth partial order matrix (adjacency matrix)

- `eps`: True error rate matrix (rows = datasets, cols = event types)

- `eps_obs`: Observed error rate matrix calculated from simulated data

- `lambdas`: Event occurrence rates (inverse of mean occurrence times)

- `T_sampling`: Sampling times for each sample (time when the sample was taken)

- `T_events`: Individual event occurrence times (N samples × n events)

- `T_sum_events`: Cumulative event times considering poset constraints (earliest possible occurrence time)

  

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

#### Please check example\panCancerCBNandMH-Sampling.R


#### Step 1: Prepare Input Data

Convert your genotype data to the required matrix format:

```r
# Example: Extract driver genes and chromosome information
# Assuming Genotype is a data frame with columns for each gene/mutation
driver_genes <- c("TP53", "KRAS", "EGFR", "MYC")  # List of driver genes
driver_chr <- c("chr17", "chr12", "chr7", "chr8")  # Corresponding chromosomes

# Convert to matrix format (samples as rows, events as columns)
mat <- as.matrix(Genotype[, c(driver_genes, driver_chr)])

# Alternatively, use simulated data from simulation_Data_Generate()
# mat <- sim_data$obs_events  # From simulated data
# Ensure the matrix contains integer values (0/1)
```


#### Step 2: Prepare Input Matrices

For multi-dataset analysis with different event types:

```r
# Get unique datasets
dataSet_unique <- unique(Genotype[,"dataSet"])

# Configure isF matrix (data filled indicator)
isF <- matrix(as.integer(c(1)), nrow = length(dataSet_unique), ncol = 2)

# Configure setD matrix (dataset ranges, 0-indexed)
setD <- as.data.frame(array(NA, dim = c(length(dataSet_unique), 2)))
for (l in 1:length(dataSet_unique)) {
  setD[l, 1] <- min(which(Genotype[,"dataSet"] == dataSet_unique[l]) - 1)
  setD[l, 2] <- max(which(Genotype[,"dataSet"] == dataSet_unique[l]) - 1)
}
setD <- as.matrix(setD)
for (l in 1:dim(setD)[2]) {
  setD[, l] <- as.integer(setD[, l])
}

# Configure eventD matrix (event type ranges)
# Event types: driver_genes (type 1) and driver_chr (type 2)
eventD <- as.data.frame(array(NA, dim = c(2, 2)))
eventD[1, 1] <- 1  # Start index of driver genes
eventD[1, 2] <- length(driver_genes)  # End index of driver genes
eventD[2, 1] <- length(driver_genes) + 1  # Start index of driver chromosomes
eventD[2, 2] <- length(driver_genes) + length(driver_chr)  # End index
eventD <- as.matrix(eventD)
for (l in 1:dim(eventD)[2]) {
  eventD[, l] <- as.integer(eventD[, l])
}

# Prepare final genotype matrix with first column as all 1s
Geno <- cbind(rep(as.integer(1), dim(mat)[1]), mat)
Geno <- as.data.frame(Geno)
for (l in 1:dim(Geno)[2]) {
  Geno[, l] <- as.integer(Geno[, l])
}
Geno <- as.matrix(Geno)

# Set poset size parameter
poset_np <- min(length(c(driver_genes, driver_chr)), 8)
```


#### Step 3: Estimate Error Rates (epsilon)

```r
# Multi-dataset epsilon estimation using estimate_Epsilon_ForMulti()
eps <- estimate_Epsilon_ForMulti(
  pat = Geno,                    # Observed genotypes matrix (rows=samples, cols=events), first column must be all 1s
  isF = isF,                     # Matrix indicating whether data in dataset[i] and eventset[j] is filled (1=True, 0=False)
  setD = setD,                   # Dataset index description matrix (rows=datasets, cols=[start, end])
  eventD = eventD,               # Eventset index description matrix (rows=eventsets, cols=[start, end])
  multi_thrds = 1,               # Number of threads for multi-dataset parallel processing
  threshold = 0.05,              # Threshold for epsilon estimate
  threshold1 = as.integer(8),    # Threshold for epsilon estimation method selection (default: 7)
  n_p = as.integer(8),           # Max number of events used for epsilon estimation
  T = 10.0,                      # Temperature of simulated annealing algorithm
  N_iter = 200L,                 # Max iteration number of simulated annealing algorithm
  lambdaS = 1.0,                 # Rate of the sampling process
  thrds = as.integer(10)         # Number of threads for parallel execution
)
```

#### Step 4: Infer Partial Order (Poset)

```r
# Find consensus poset using voting with find_Poset_ForVote()
poset <- find_Poset_ForVote(
  pat = Geno,                    # Observed genotypes matrix
  isF = isF,                     # Matrix indicating whether data is filled
  eps = eps,                     # Error rate matrix from Step 3
  setD = setD,                   # Dataset index description matrix
  eventD = eventD,               # Eventset index description matrix
  Fine_Tune_Num = as.integer(2L),# Fine tune number for more accurate poset (best=2)
  vote_size = 5L,                # Number of votes to find consensus poset
  vote_threshold = 0.5,          # Threshold in the vote (0-1)
  vote_thrds = 1L,               # Number of threads in the vote
  threshold = 0.0001,            # Threshold for determining partial order relationships (filter false positives)
  threshold2 = 0.01,             # Threshold considering tolerance (should be >= threshold)
  n_p = as.integer(poset_np),    # Number of events during fine-tune
  is_update_eps = FALSE,         # Whether update eps in fine tune (experience shows FALSE is better)
  T = 10,                        # Temperature of simulated annealing algorithm
  N_iter = 200L,                 # Max iteration number of simulated annealing algorithm
  lambdaS = 1.0,                 # Rate of the sampling process
  thrds = as.integer(20)         # Number of threads for parallel execution
)
```

#### Step 5: Estimate Rate Parameters (lambda)

```r
# Estimate lambda using EM algorithm with estimate_Lambda()
fit <- estimate_Lambda(
  pat = Geno,                    # Observed genotypes matrix
  poset = poset,                 # Partial order adjacency matrix from Step 4
  isF = isF,                     # Matrix indicating whether data is filled
  eps = eps,                     # Error rate matrix
  setD = setD,                   # Dataset index description matrix
  eventD = eventD,               # Eventset index description matrix
  lambdaS = 1.0,                 # Rate of the sampling process
  L = 100L,                      # Number of samples drawn from proposal in E-step
  sampling = 'add-remove',       # Sampling scheme: "forward", "add-remove", "backward", "bernoulli", or "pool"
  maxIter = 500L,                # Maximum number of EM iterations
  updateStepSize = 20L,          # EM steps after which L is doubled if convergence not reached
  tol = 0.001,                   # Convergence tolerance for error rate and rate parameters
  maxLambda = 1e6,               # Upper bound on the value of rate parameters
  neighborhoodDist = 1L,         # Hamming distance for "backward" sampling
  is_update_eps = FALSE,         # Whether update eps during lambda estimation (experience shows FALSE is better)
  thrds = as.integer(20)         # Number of threads for parallel execution
)

# Store additional results
fit$lambdaS <- 1.0
fit$poset <- poset
```

## 5. Extended Analysis (Project-Specific)

After obtaining the CBN results (poset, epsilon, lambda) from the basic workflow, you can perform additional project-specific analyses using the scripts in the `example/` directory. These analyses are tailored to specific research questions and require additional data sources.

### 5.1 Absolute Time Calibration with MH Sampling

**Script**: `example/panCancerCBNandMH-Sampling.R`

### 5.2 Immune Fitness Inference

**Script**: `example/panCancerImmuneFitnessInference.R`

### 5.3 Single-Cell RNA Analysis

**Script**: `example/scRNA_Analysis.R`


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

| Function                                       | Description                                |
| ---------------------------------------------- | ------------------------------------------ |
| `estimate_Epsilon()`                           | Estimate error rates for hidden CBN model  |
| `estimate_Epsilon_based_on_Poset_and_Lambda()` | Epsilon estimation with known poset/lambda |
| `estimate_Epsilon_ForMulti()`                  | Multi-dataset epsilon estimation           |
| `find_Poset()`                                 | Infer partial order from genotype data     |
| `find_Poset_ForVote()`                         | Consensus poset via voting                 |
| `estimate_Lambda()`                            | Estimate event rate parameters             |
| `estimate_Lambda_ForMulti()`                   | Multi-dataset lambda estimation            |

### Utility Functions

| Function                                      | Description                       |
| --------------------------------------------- | --------------------------------- |
| `topological_Sort()`                          | Perform topological sort on poset |
| `random_Poset()`                              | Generate random poset structure   |
| `rateTimeHelp()`                              | Convert lambda to time bounds     |
| `find_most_Compatible_Genotype_by_Flipping()` | Find compatible genotypes         |
| `donor_Pair_Genotype_Filter()`                | Filter donor pair genotypes       |
| `GammaCluster()`                              | Gaussian Mixture Model clustering |

### Data Generation Functions

| Function                          | Description                      |
| --------------------------------- | -------------------------------- |
| `simulation_Data_Generate()`      | Generate simulated genotype data |
| `simulation_Time_Data_Generate()` | Generate simulated timing data   |

---




## Support

For questions and support, contact: Hu Changbao <1437894182@qq.com>

## License

GPL (≥ 2)
