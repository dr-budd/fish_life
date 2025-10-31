# :tropical_fish::dna: Fish lifepsan prediction from genomic data

## :page_facing_up: Associated Publication

[![Journal](https://img.shields.io/badge/Published_in-Molecular_Ecology_Resources-0072B8.svg)](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13774)

> **Budd, A.M.**, Mayne, B., Berry, O. & Jarman, S. (2025). Fish species lifespan prediction from promoter cytosine-phosphate-guanine density. *Molecular Ecology Resources*, **25**, e13774. [https://doi.org/10.1111/1755-0998.13774](https://doi.org/10.1111/1755-0998.13774)

---

## :file_folder: Structure

The code is separated into the following scripts that should be run in order:

- `01_wget_genomes.sh`
- `02_compare_checksums.sh`
- `03_unzip_genomes.slurm`
- `04_create_blast_dbs.slurm`
- `05_query_blast_dbs.slurm`
- `06_calculate_cpg_content.R`
- `06_calculate_cpg_content_fc.R`
- `07_final_data_set.Rmd`
- `08_elastic_net_nested_cv_oe.R`
- `09_elastic_net_nested_cv_results_bagged.Rmd`
- `10_lifespan_model_bagged.Rmd`
- `11_promoter_functional_anal_gprofiler.R`
- `12_tree_figure.R`

Please note, any bash or slurm scripts (`.sh` or `.slurm` extension) were written for and run on a HPC facility that uses a SLURM batch-queue system. This means that many of the slurm scripts specify core allocation, run times and memory usage allocation that may need to be adapted for different platforms.

## :chart_with_upwards_trend: Data

The raw genomic data which can be downloaded directly from https://www.ncbi.nlm.nih.gov/genome/ using the provided accession numbers. Additional files are in the folder `dataFiles`.

## :woman_technologist: Author
Alyssa Budd (alyssa.budd@csiro.au)

## :copyright: License
[CSIRO Open Source Software Licence Agreement (variation of the BSD / MIT License)](LICENSE.txt)
