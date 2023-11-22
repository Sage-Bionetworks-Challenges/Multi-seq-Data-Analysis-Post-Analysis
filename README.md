# Multi-seq-Data-Analysis-Post-Analysis

Post-hoc analysis for [scRNA-seq and scATAC-seq Data Analysis DREAM Challenge](https://www.synapse.org/#!Synapse:syn26720920/wiki/615338)

## Installation

1.  Clone the repo

        git clone https://github.com/Sage-Bionetworks-Challenges/Multi-seq-Data-Analysis-Post-Analysis
        cd Multi-seq-Data-Analysis-Post-Analysis

2.  Create a [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) environment using python 3.9:

        conda create --name synapse python=3.9 -y
        conda activate synapse

3.  Install Python dependencies

        python -m pip install challengeutils==4.2.0

    check if `synapseclient` and `challengeutils` are installed via:

        synapse --version
        challengeutils -v

4.  Install R dependencies

        R -e 'source("install.R")'

    > **Note:** <br>
    > The task 2 analysis uses `bedr` package that has two requisitions - [bedpos](https://anaconda.org/bioconda/bedops) and [tabix](https://anaconda.org/bioconda/tabix) needed to be installed as well.

5.  Set up Synapse credentials via CLI, or manually store the credentials to `~/.synapseConfig` - see details [here](https://help.synapse.org/docs/Client-Configuration.1985446156.html)
    synapse login --rememberMe

## Usage

Download all final submission results and each individual test case's scores to `data/` folder:

    Rscript submission/get_submissions.R

- `final_submissions_{task}.rds`: Esseential information of final submission, e.g submission id, team, ranks
- `final_scores_{task}.rds`: All test case scores from final submissions, consists of test case name, scores of primary and secondary metrics

Download output files (imputed gene expression / called peaks) of all final submissions to `data/model_output/`

    # replace {task} with 'task1' or 'task2'
    Rscript submission/get_predictions_{task}.R

> **Warning**
> For Task 1, the output (imputation) of each submission has large size ~30G. Please be aware of the available disk space.

Report statistics about submissions
    
    Rscript -e 'rmarkdown::render("stats/get_submission_stats.rmd")'
