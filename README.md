# Multi-seq-Data-Analysis-Post-Analysis

Post-hoc analysis for scRNA-seq and scATAC-seq Data Analysis DREAM Challenge

## Set up environment

    conda create --name synapse python=3.9 -y
    conda activate synapse

## Install Python dependecies

    python -m pip install challengeutils==4.2.0

check if `synapseclient` and `challengeutils` are installed via:

    synapse --version
    challengeutils -v 

## Install R dependencies

    R -e 'source("install.R")'

> **Note** install \`bedr\`'s requisitions [bedpos](https://anaconda.org/bioconda/bedops) and [tabix](https://anaconda.org/bioconda/tabix) as well for Task 2
