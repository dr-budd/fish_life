## SET UP ----

## set working directly to query results directory
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/dataFiles/queryResultsFc"))

## open file to write output
sink("../../LogFiles/06_calculate_cpg_content_fc_log.txt")

## load libraries
library(stringr)
library(ape)
library(tidyverse)

## create a list of blast hit files to import into R
data_files_all <- list.files(pattern = "top_hits")

## remove any files of size zero from file list
data_files_empty <- NULL
for (i in data_files_all) {
  if (file.size(i) == 0)
    data_files_empty <- c(data_files_empty, i)
}

data_files <- NULL
for (i in data_files_all) {
  if (!file.size(i) == 0)
    data_files <- c(data_files, i)
}

## create vector of assembly names
assembly_names <- gsub("_genomic_top_hits", "", data_files)

## CALCULATE CpG DENSITY ----

## for first file in dir ----

## specify the file number
file_num <- 1

## save the assembly name to a variable (the name of the file minus the extension)
temp_assembly <- paste(gsub("_genomic_top_hits", "", data_files[file_num]))

## print the file number and assembly name
print(paste("file number:", file_num, "of", length(data_files)))

## import and formate data, then calculate CG density
temp_data <- read.table(file = data_files[file_num]) %>%
  ## add column names (from blast query call)
  ## NB: this is necessary but nice if you want to eyeball the data
  rename(query_seqid = V1, subject_seqid = V2, subject_seq_len = V3, 
         query_start = V4, query_end = V5, ## start and end of alignment in query
         subject_start = V6, subject_end = V7, ## start and end of alignment in subject
         perc_ident = V8, num_ident = V9, ## percent and number of identical matches
         subject_seq = V10) %>% ## aligned part of subject sequence 
  ## remove sequence gap indicators from alignment ("-")
  mutate(subject_seq = gsub("-", "", subject_seq), 
         ## remove sequence gap indicators from alignment ("-")
         length = nchar(subject_seq), 
         ## count the number of "CG" occurrences
         cg_count = str_count(subject_seq, "CG"), 
         ## count the number of Cs
         c_count = str_count(subject_seq, "C"),
         ## count the number of Gs
         g_count = str_count(subject_seq, "G"),
         ## calculate the CpG density
         cg_dens = cg_count/length,
         ## calculate the CpG O/E 
         cg_oe = cg_count/((c_count*g_count)/length)) %>%
  ## rename the promoter ID column
  rename(., promoter_id = query_seqid)

## remove duplicate hits (retain only 'top' hit)
temp_data <- temp_data[!duplicated(temp_data[,1]), ]

## edit and rename the OE data frame 
temp_oe_data <- temp_data %>%
  ## retain only the promoter id and CpG densities
  select(promoter_id, cg_oe) %>%
  ## rename the cg_oe column with the assembly name
  rename(., "{temp_assembly}" := cg_oe) %>%
  ## reorder by promoter id
  arrange(., promoter_id)

## edit and rename the density data frame 
temp_dens_data <- temp_data %>%
  ## retain only the promoter id and CpG densities
  select(promoter_id, cg_dens) %>%
  ## rename the cg_oe column with the assembly name
  rename(., "{temp_assembly}" := cg_dens) %>%
  ## reorder by promoter id
  arrange(., promoter_id)

## edit and rename the length data frame 
temp_len_data <- temp_data %>%
  ## retain only the promoter id and hit lengths
  select(promoter_id, length) %>%
  ## rename the cg_oe column with the assembly name
  rename(., "{temp_assembly}" := length) %>%
  ## reorder by promoter id (optional here)
  arrange(., promoter_id)

## edit and rename the percent id data frame 
temp_pi_data <- temp_data %>%
  ## retain only the promoter id and percent idents
  select(promoter_id, perc_ident) %>%
  ## rename the cg_oe column with the assembly name
  rename(., "{temp_assembly}" := perc_ident) %>%
  ## reorder by promoter id (optional here)
  arrange(., promoter_id)

## for remaining files in dir ----

## loop through each file 
for(file_num in 2:length(data_files)){
  
  ## save the assembly name to a variable (the name of the file minus the extension)
  temp_assembly <- paste(gsub("_genomic_top_hits", "", data_files[file_num]))
  
  ## print the file number and assembly name
  print(paste("file number:", file_num, "of", length(data_files)))
  
  ## import and format data, then calculate CG density
  temp_data <- read.table(file = data_files[file_num]) %>%
    ## add column names (from blast query call)
    ## NB: this is not necessary but nice if you want to eyeball the data
    rename(query_seqid = V1, subject_seqid = V2, subject_seq_len = V3, 
           query_start = V4, query_end = V5, ## start and end of alignment in query
           subject_start = V6, subject_end = V7, ## start and end of alignment in subject
           perc_ident = V8, num_ident = V9, ## percent and number of identical matches
           subject_seq = V10) %>% ## aligned part of subject sequence 
    ## remove sequence gap indicators from alignment ("-")
    mutate(subject_seq = gsub("-", "", subject_seq), 
           ## remove sequence gap indicators from alignment ("-")
           length = nchar(subject_seq), 
           ## count the number of "CG" occurrences
           cg_count = str_count(subject_seq, "CG"), 
           ## count the number of Cs
           c_count = str_count(subject_seq, "C"),
           ## count the number of Gs
           g_count = str_count(subject_seq, "G"),
           ## calculate the CpG density
           cg_dens = cg_count/length,
           ## calculate the CpG O/E 
           cg_oe = cg_count/((c_count*g_count)/length)) %>%
    ## rename the promoter ID column
    rename(., promoter_id = query_seqid)
  
  ## remove duplicate hits (retain only 'top' hit)
  temp_data <- temp_data[!duplicated(temp_data[,1]), ]
  
  ## edit and rename the OE data frame 
  temp_oe_data <- temp_data %>%
    ## retain only the promoter id and CpG densities
    select(promoter_id, cg_oe) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := cg_oe) %>%
    ## reorder by promoter id
    arrange(., promoter_id) %>%
    ## join temp data to existing data (and overwrite it)
    full_join(temp_oe_data, ., by = "promoter_id") 
  
  ## edit and rename the density data frame 
  temp_dens_data <- temp_data %>%
    ## retain only the promoter id and CpG densities
    select(promoter_id, cg_dens) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := cg_dens) %>%
    ## reorder by promoter id
    arrange(., promoter_id) %>%
    ## join temp data to existing data (and overwrite it)
    full_join(temp_dens_data, ., by = "promoter_id") 
  
  ## join the data frames (hit length)
  temp_len_data <- temp_data %>%
    ## retain only the promoter id and hit lengths
    select(promoter_id, length) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := length) %>%
    ## reorder by promoter id (optional here)
    arrange(., promoter_id) %>%
    ## join temp data to existing data (and overwrite it)
    full_join(temp_len_data, ., by = "promoter_id") 
  
  ## edit and rename the percent id data frame 
  temp_pi_data <- temp_data %>%
    ## retain only the promoter id and percent idents
    select(promoter_id, perc_ident) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := perc_ident) %>%
    ## reorder by promoter id (optional here)
    arrange(., promoter_id) %>%
    ## join
    full_join(temp_pi_data, ., by="promoter_id")
  
}

## format cpg oe data  
cpg_oe_data_full <- temp_oe_data %>%
  ## order the dataframe by promoter id 
  arrange(., promoter_id) %>%
  ## transpose CpG oe data for elastic net
  ## this is round about but prevents num to chr 
  pivot_longer(-promoter_id) %>%
  pivot_wider(names_from = promoter_id, values_from = value) %>%
  ## move the assembly names to rownames
  column_to_rownames(., "name") %>%
  ## convert to matrix
  as.matrix(.)

## export cpg oe file (all fish)
write.csv(cpg_oe_data_full, "../06.00_cpg_oe_data_full_fc.csv", row.names = TRUE)

## format cpg density data
cpg_dens_data_full <- temp_dens_data %>%
  ## order the dataframe by promoter id
  arrange(., promoter_id) %>%
  ## transpose CpG density data for elastic net
  ## this is round about but prevents num to chr
  pivot_longer(-promoter_id) %>%
  pivot_wider(names_from = promoter_id, values_from = value) %>%
  ## move the assembly names to rownames
  column_to_rownames(., "name") %>%
  ## convert to matrix
  as.matrix(.)

## export cpg density file (all fish)
write.csv(cpg_dens_data_full, "../06.00_cpg_dens_data_full_fc.csv", row.names = TRUE)

## format hit length data
hit_len_data_full <- temp_len_data %>%
  ## order the dataframe by promoter id 
  arrange(., promoter_id) %>%
  ## transpose CpG density data for elastic net
  ## this is round about but prevents num to chr 
  pivot_longer(-promoter_id) %>%
  pivot_wider(names_from = promoter_id, values_from = value) %>%
  ## move the assembly names to rownames
  column_to_rownames(., "name") %>%
  ## convert to matrix
  as.matrix(.)

## export hit length file (all fish)
write.csv(hit_len_data_full, "../06.00_hit_len_data_full_fc.csv", row.names = TRUE)

## format percent identity data
hit_pi_data_full <- temp_pi_data %>%
  ## order the dataframe by promoter id 
  arrange(., promoter_id) %>%
  ## transpose CpG density data for elastic net
  ## this is round about but prevents num to chr 
  pivot_longer(-promoter_id) %>%
  pivot_wider(names_from = promoter_id, values_from = value) %>%
  ## move the assembly names to rownames
  column_to_rownames(., "name") %>%
  ## convert to matrix
  as.matrix(.)

## export cpg density file (all fish)
write.csv(hit_pi_data_full, "../06.00_hit_pi_data_full_fc.csv", row.names = TRUE)

## CREATE METADATA ----

## calculate blast hit info ----

## calculate mean hit length per species assembly
nz_hit_len_all_fish <- temp_len_data %>%
  ## order the data frame by promoter id 
  arrange(., promoter_id) %>%
  ## format to long and get summary data
  pivot_longer(-promoter_id, names_to = "assembly_id") %>%
  group_by(assembly_id) %>%
  summarise(mean_nz_hit_len = mean(value, na.rm = TRUE))

## calculate mean percent identity per species assembly
nz_pi_all_fish <- temp_pi_data %>%
  ## order the dataframe by promoter id 
  arrange(., promoter_id) %>%
  ## format to long and get summary data
  pivot_longer(-promoter_id, names_to = "assembly_id") %>%
  group_by(assembly_id) %>%
  summarise(mean_nz_perc_ident = mean(value, na.rm = TRUE))

## calculate mean hit length per species assembly
hit_len_all_fish <- temp_len_data %>%
  ## replace NAs with zeros
  replace(is.na(.), 0) %>%
  ## order the dataframe by promoter id
  arrange(., promoter_id) %>%
  ## format to long and get summary data
  pivot_longer(-promoter_id, names_to = "assembly_id") %>%
  group_by(assembly_id) %>%
  summarise(mean_hit_len = mean(value, na.rm = TRUE))

## calculate mean hit length per species assembly
pi_all_fish <- temp_pi_data %>%
  ## replace NAs with zeros
  replace(is.na(.), 0) %>%
  ## order the dataframe by promoter id
  arrange(., promoter_id) %>%
  ## format to long and get summary data
  pivot_longer(-promoter_id, names_to = "assembly_id") %>%
  group_by(assembly_id) %>%
  summarise(mean_perc_ident = mean(value, na.rm = TRUE))

## calculate number of hits per species assembly 
hit_num_all_fish <- as.data.frame(cpg_oe_data_full) %>% 
  mutate(num_hits = rowSums(!is.na(.))) %>%
  select(num_hits) %>%
  rownames_to_column(., var = "assembly_id")

## import metadata sum and add blast hit info to metadata 
metadata_sum <- read.csv("../00.11_metadata_sum_fc.csv") %>%
  ## add non-zero hit length
  left_join(., nz_hit_len_all_fish, by = "assembly_id") %>%
  ## add hit length
  left_join(., hit_len_all_fish, by = "assembly_id") %>%
  ## add num hits
  left_join(., hit_num_all_fish, by = "assembly_id") %>%
  ## add perc id
  left_join(., pi_all_fish, by = "assembly_id") %>%
  ## add non-zero perc id
  left_join(., nz_pi_all_fish, by = "assembly_id") 

## SUBSET METADATA ----

## remove species assemblies in metadata but not cg oe data
no_cgoe_fish <- setdiff(metadata_sum$assembly_id, rownames(cpg_oe_data_full))
paste0("removing assembly with no cpg oe data from the metadata: ", no_cgoe_fish)
metadata_sum <- metadata_sum %>%
  filter(., !assembly_id %in% no_cgoe_fish)

## remove species assemblies in the cg oe and dens data but not the metadata
no_meta_fish <- setdiff(rownames(cpg_oe_data_full), metadata_sum$assembly_id)
paste0("removing assembly with no lifespan info from the cpg oe and cg dens data: ", no_meta_fish)
cpg_oe_data <- cpg_oe_data_full[!rownames(cpg_oe_data_full) %in% no_meta_fish, ]
cpg_dens_data <- cpg_dens_data_full[!rownames(cpg_dens_data_full) %in% no_meta_fish, ]

## add mean hit info to metadata
metadata_sum <- metadata_sum %>%
  ## add mean cpg oe
  left_join(., (as.data.frame(cpg_oe_data) %>%
                  ## replace NAs with zeros
                  replace(is.na(.), 0)) %>%
              mutate(mean_cpg_oe = rowMeans(.)) %>%
              select(mean_cpg_oe) %>%
              rownames_to_column(., var = "assembly_id"),
            by = "assembly_id") %>%
  ## add mean non-zero cpg oe
  left_join(., (as.data.frame(cpg_oe_data)) %>%
              mutate(mean_nz_cpg_oe = rowMeans(., na.rm=TRUE)) %>%
              select(mean_nz_cpg_oe) %>%
              rownames_to_column(., var = "assembly_id"),
            by = "assembly_id") %>%
  ## add mean cpg dens
  left_join(., (as.data.frame(cpg_dens_data) %>%
                  ## replace NAs with zeros
                  replace(is.na(.), 0)) %>%
              mutate(mean_cpg_dens = rowMeans(.)) %>%
              select(mean_cpg_dens) %>%
              rownames_to_column(., var = "assembly_id"),
            by = "assembly_id") %>%
  ## add mean non-zero cpg dens
  left_join(., (as.data.frame(cpg_dens_data)) %>%
              mutate(mean_nz_cpg_dens = rowMeans(., na.rm=TRUE)) %>%
              select(mean_nz_cpg_dens) %>%
              rownames_to_column(., var = "assembly_id"),
            by = "assembly_id") %>%
  ## arrange by assembly id
  arrange(., assembly_id) 

## quick check
ifelse(all(rownames(cpg_oe_data) == metadata_sum$assembly_id) == "TRUE",
       "assembly ids are identical - you are good to go",
       paste0("STOP!!! something is wrong - your assembly ids are not identical in your cpg oe and your metadata (and they need to be)"))

## EXPORT FILES ----

## export cpg oe file
write.csv(cpg_oe_data, "../06.00_cpg_oe_data_fc.csv", row.names = TRUE)

## export cpg density file
write.csv(cpg_dens_data, "../06.00_cpg_dens_data_fc.csv", row.names = TRUE)

## export metadata
write.csv(metadata_sum, "../06.00_metadata_sum_fc.csv", row.names = FALSE)

## print session info
sessionInfo()

## write output to open file
sink()

## END SCRIPT
##
