HIV Transmission Dynamics in Greater Mexico City

This repository brings together data and analysis workflows supporting research on HIV transmission dynamics and spatial mixing in Greater Mexico City, Mexico. It combines phylogenetic and phylogeographic analyses, incidence-informed sequence subsampling, and mobility-network analyses for Mexico City (CDMX) and the State of Mexico (MEX).

Repository contents

•
Processing_phylogenetic_data/ — R workflows, phylogenetic trees, BEAST analysis templates, Mexican geographic files, population data, and derived inputs used to investigate introductions, clusters, spatial dispersal, and mobility-related patterns.

•
Subsample_Incidence/ — Python scripts and supporting files for incidence-informed HIV sequence subsampling and year-by-location sequence counts, including analyses of potential geographic sampling bias.

•
mobility_networks/ — Python scripts for constructing daily AGEB-level origin–destination networks, summarizing directional mobility between CDMX and the State of Mexico, mapping daily flows, and describing daily device and record counts.

Each folder contains its own scripts and, where available, detailed documentation. Some workflows rely on external source data or project-specific file paths; review the folder-level README and update input/output paths before running them. The files in this repository include derived datasets and analysis materials, but may not include every raw dataset required to reproduce all results.

