# :tropical_fish::dna: Fish lifepsan prediction from genomic data

Analysis pipeline for:

>Budd, A.M., Mayne, B., Berry, O. and Jarman, S., 2023. Fish species lifespan prediction from promoter cytosine‐phosphate‐guanine density. Molecular Ecology Resources.

:memo: https://doi.org/10.1111/1755-0998.13774

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

The raw genomic data which can be downloaded directly from https://www.ncbi.nlm.nih.gov/genome/ using the provided accession number. Additional files are in the folder `dataFiles`.

## :woman_technologist: Author
Alyssa Budd (alyssa.budd@csiro.au)

## :bouquet: Acknowledgements
The scripts have been adapted from a previous model used to predict fish lifespan https://github.com/dr-budd/fish_life. All scripts in `lifepi` were by written by Suk Yee Yong (sukyee.yong@csiro.au).

## :copyright: License
[CSIRO Open Source Software Licence Agreement (variation of the BSD / MIT License)](LICENSE.txt)
