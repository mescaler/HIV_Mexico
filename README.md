# HIV Transmission Dynamics in Greater Mexico City

This repository contains the data, scripts, and analytical workflows supporting our research on **HIV transmission dynamics and spatial mixing in Greater Mexico City, Mexico**.

The study integrates phylogenetic and phylogeographic analyses, incidence-informed sequence subsampling, and human mobility networks to investigate HIV transmission patterns across Mexico City (CDMX) and the State of Mexico (MEX).

## Repository Structure

The repository is organised into three main directories:

### 1. `Processing_phylogenetic_data/`
**Phylogenetic and phylogeographic analyses**

Contains R workflows and supporting materials for investigating HIV evolutionary and spatial dynamics.

Includes:
- BEAST xmls
- Scripts for processing phylogenetic and phylogeographic results
- Derived outputs for investigating viral introductions, genetic clusters, spatial dispersal, and mobility-related patterns

### 2. `Subsample_Incidence/`
**Incidence-informed subsampling and sampling bias assessment**

Contains Python scripts for evaluating the potential influence of geographic sampling imbalances on phylogeographic inference.
- Incidence-informed HIV sequence subsampling
- Supporting files for assessing geographic sampling bias

### 3. `mobility_networks/`
**Human mobility network analysis**

Contains Python scripts for processing and analysing human mobility data across Greater Mexico City.

Includes:
- Construction of daily origin-destination networks at the AGEB level
- Quantification of directional mobility between CDMX and MEX
- Mapping and visualisation of daily mobility flows
- Summaries of daily device and mobility record counts

**Before running the analyses:**
1. Review the documentation within each directory.
2. Update project-specific input and output file paths.
3. Ensure that the required external datasets and software dependencies are available.

Some workflows depend on external source data that are not distributed within this repository (see SI of original manuscript). 
