# FMS_B_cell_tolerance

Code accompanying the manuscript **“Defective peripheral B cell tolerance leads to dysregulated B cell responses in Fibromyalgia Syndrome”**.

This repository contains R scripts used to analyse B cell receptor (BCR) repertoire sequencing and flow-cytometry data from individuals with fibromyalgia syndrome (FMS) and healthy controls, and to generate the analyses and figures reported in the manuscript. The repository includes BCR repertoire analysis, flow-cytometry stimulation/proliferation analysis, statistical testing, and figure generation. The repository is released under the GPL-3.0 licence.

## 1. System requirements

### Operating system

The analyses were developed and run using **R 4.2.3** for the flow-cytometry analyses. The code is written in R and should in principle run on Windows, macOS, or Linux, provided the required R packages and input files are available.

**Tested software environment**

* R 4.2.3
* RStudio: not required
* Git: recommended for obtaining the repository

The repository does not currently provide a fully locked package environment (e.g. `renv.lock`). Exact package versions should therefore be recorded from the analysis environment used to generate the manuscript results before attempting exact computational reproduction.

### R packages

The analyses use standard R/statistical and plotting packages, including packages for:

* flow-cytometry processing and transformation
* statistical modelling
* data manipulation
* plotting
* heatmaps
* proliferation modelling
* PCA and multivariate analysis

The flow-cytometry analysis specifically uses `flowCore`, `ggplot2`, `MASS`, and `nor1mix`. The manuscript also describes use of `pheatmap` and standard/base R functions. Mixed-effects analyses described for the single-cell analysis use `lme4` and `lmerTest`.

The BCR analysis uses custom R scripts together with previously processed BCR repertoire input tables.

### Hardware

No specialised hardware is required. A standard desktop or laptop computer with at least **8 GB RAM** should be sufficient for most analyses. Larger memory may be beneficial for processing large BCR repertoire tables or flow-cytometry datasets.

No GPU is required.

## 2. Installation guide

### Obtain the repository

Clone the repository using Git:

```bash
git clone https://github.com/rbr1/FMS_B_cell_tolerance.git
cd FMS_B_cell_tolerance
```

Alternatively, download the repository as a ZIP archive from GitHub and extract it locally.

### Install R

Install R 4.2.3 or a later compatible R release.

Install the required R packages, for example:

```r
install.packages(c(
  "ggplot2",
  "MASS",
  "pheatmap",
  "lme4",
  "lmerTest",
  "nor1mix"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
```

Additional packages used by individual scripts may need to be installed according to the script being run.

### Installation time

On a typical modern desktop computer with internet access, installation of R and the required packages should generally take **5–20 minutes**, depending on the operating system and whether packages need to be compiled from source.

## 3. Demo

The repository is primarily organised for reproduction of the analyses performed for the manuscript rather than as a standalone command-line application. The scripts currently expect input data and directory structures corresponding to the original analysis environment.

### BCR analysis

The main BCR analysis scripts are:

```text
BCR_Print raw_values.R
BCR_Raw_input_value_Overall analysis_MANOVA.R
```

These scripts analyse processed BCR repertoire tables and generate repertoire-level measurements, statistical comparisons, and figures.

Example:

```r
source("BCR_Print raw_values.R")
source("BCR_Raw_input_value_Overall analysis_MANOVA.R")
```

Before running, update the input and output directory variables in the scripts to point to the location of the corresponding BCR input files and desired output directory.

The analysis includes:

* immunoglobulin isotype usage
* IGHV and IGHJ gene usage
* clonality measures
* somatic hypermutation
* CDR3 characteristics
* repertoire comparisons between FMS and healthy controls
* multivariate statistical analyses
* PCA and associated plots

### Flow-cytometry analysis

The main flow-cytometry scripts are:

```text
FACS_stimulation_analysis.R
FACS_stimulation_proliferation plot.R
```

These analyse compensated/transformed flow-cytometry data from unstimulated and stimulated PBMC samples and generate analyses of:

* B cell subset composition
* CD21/CD22/CD24 phenotypes
* proliferation indices
* division indices
* stimulation responses
* associated figures

Input flow-cytometry files are expected in the directory structure specified within the scripts and accompanying `File_information.txt`.

### Expected output

Outputs are primarily:

* PDF figures
* statistical summaries
* intermediate analysis tables
* PCA and repertoire plots
* proliferation and flow-cytometry visualisations

The exact filenames and locations are determined by the output-directory settings in each script.

### Expected runtime

Runtime depends strongly on the input dataset and storage location.

As a guide:

* individual figure-generation scripts: **seconds to a few minutes**
* BCR repertoire processing and statistical analysis: **minutes to tens of minutes**
* large flow-cytometry analyses: **minutes to longer**, depending on number of FCS files and available memory

Because the repository is not currently packaged as a single automated workflow, exact end-to-end runtime has not been benchmarked.

## 4. Instructions for use

### Using the code with new data

The repository is intended primarily for reproducing and extending the analyses reported in the manuscript.

For BCR repertoire data:

1. Prepare processed BCR repertoire tables in the format expected by the BCR scripts.
2. Update the input/output directory variables in:

   * `BCR_Print raw_values.R`
   * `BCR_Raw_input_value_Overall analysis_MANOVA.R`
3. Run the scripts sequentially.
4. Inspect the generated statistical summaries and PDF figures.

The scripts expect BCR annotations including sample identifiers, immunoglobulin isotype, IGHV/IGHJ gene assignments, sequence counts, and clone-related measurements.

For flow-cytometry data:

1. Organise FCS files according to the directory structure expected by the analysis script.
2. Update the data and output paths.
3. Ensure that the relevant FMO controls and threshold files are available.
4. Run:

   * `FACS_stimulation_analysis.R`
   * `FACS_stimulation_proliferation plot.R`
5. Review the generated subset, proliferation, and statistical outputs.

Threshold and gate information used by the analysis is provided in:

```text
MFI_gates.txt
MFI_thresholds.txt
Proliferation_thresholds.txt
```

`File_information.txt` documents the file paths and input files used in the original flow-cytometry analysis.

## Reproduction of manuscript analyses

The repository contains the analysis code used to generate the manuscript results and figures. The manuscript cohort included BCR repertoire data from **15 individuals with FMS and 19 healthy controls**, with flow-cytometry experiments performed in a subset of **14 FMS and 14 healthy controls**. The BCR data were generated from PBMC RNA using UMI-containing, isotype-specific primers and Illumina paired-end sequencing. Flow-cytometry data were generated from PBMCs analysed at baseline and following approximately 5 days of stimulation.

For exact reproduction, users require:

1. the analysis code in this repository;
2. the corresponding processed input data;
3. the original raw sequencing/flow-cytometry data where applicable;
4. the same or documented compatible versions of R and R packages.

BCR sequencing data and associated metadata are being deposited separately in the European Genome-phenome Archive (EGA).

## Repository structure

```text
FMS_B_cell_tolerance/
├── BCR_Print raw_values.R
├── BCR_Raw_input_value_Overall analysis_MANOVA.R
├── FACS_stimulation_analysis.R
├── FACS_stimulation_proliferation plot.R
├── File_information.txt
├── MFI_gates.txt
├── MFI_thresholds.txt
├── Proliferation_thresholds.txt
├── README.md
└── LICENSE
```

## Data availability

The repository contains analysis code and supporting configuration/threshold files. Raw human sequencing and flow-cytometry data are not included in this repository because they contain controlled-access human participant data.

BCR sequencing data will be deposited in the **European Genome-phenome Archive (EGA)**. Access to controlled datasets is subject to the relevant EGA access procedures.
