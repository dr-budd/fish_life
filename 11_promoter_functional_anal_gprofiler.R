## data_set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load libraries
library(tidyverse)
library(gprofiler2)
library(clusterProfiler)
library(enrichplot)
library(DOSE) ## needed to convert to enrichResult object
library(viridis)
library(ggpubr)

## import data
promoters <- read.table("DataFiles/epd_danio_rerio_names.txt") %>%
  rename(promoter_id = V1, name = V2) %>%
  mutate(name = gsub("_1", "", name))

## import coefficient for predictive promoters
predictive_promoters_avg <- read.csv("DataFiles/10.00_predictive_promoters_avg.csv")

## edit df
promoter_coefficients <- predictive_promoters_avg %>%
  ## remove abs
  select(-abs_mean_coef) %>%
  ## rename
  rename(coefficient = mean_coef) %>%
  ## add gene names 
  left_join(., promoters) %>%
  ## add pos and neg category
  mutate(cor_type = ifelse(coefficient < 0, 
                       "negative", 
                       "positive"))

## create vector of genes in the best model 
model_genes <- promoter_coefficients$name

## split into negative and positive 
positive_promoters = promoter_coefficients %>%
  filter(cor_type == "positive") %>%
  select(name) %>%
  unique(.) %>%
  pull(.)

negative_promoters = promoter_coefficients %>%
  filter(cor_type == "negative") %>%
  select(name) %>%
  unique(.) %>%
  pull(.)

## create vector of genes in all 10 models 
model_genes_positive <- intersect(model_genes, positive_promoters)
model_genes_negative <- intersect(model_genes, negative_promoters)

## POSITIVE ----

## enrichment analysis using gene names
gost_results <- gost(query = model_genes_positive, 
                     organism = "drerio", 
                     multi_query = FALSE, 
                     evcodes = TRUE)

# modify the g:Profiler data frame
gost_mod_pos = gost_results$result[,c("query", "source", "term_id",
                                  "term_name", "p_value", "query_size", 
                                  "intersection_size", "term_size", 
                                  "effective_domain_size", "intersection")]

names(gost_mod_pos) = c("Cluster", "Category", "ID", "Description", "p.adjust",
                    "query_size", "Count", "term_size", "effective_domain_size",
                    "geneID")

gost_mod_pos$geneID = gsub(",", "/", gost_mod_pos$geneID)

## define as enrichResult object
gost_mod_pos_enrich  = new("enrichResult", result = gost_mod_pos)

## create plot data
plot_data <- gost_mod_pos_enrich@result %>%
  rownames_to_column(., var = "go_term") %>%
  mutate(neg_log10_padjust = -log10(p.adjust), 
         Category = gsub("GO:BP", "Biological process", Category) %>%
           gsub("GO:CC", "Cellular component", .) %>%
           gsub("GO:MF", "Molecular function", .)) %>%
  ## arrange by pvalue
  arrange(., desc(neg_log10_padjust)) 

## make category and ordered factor
plot_data$Description = factor(plot_data$Description, levels = plot_data$Description)

## plot as per bens suggestion
positive_plot <- ggplot(plot_data, aes(x=neg_log10_padjust, 
                                       fill = Category,
                                       y=Description))+
  geom_bar(stat="identity")+
  geom_text(aes(label = paste0("n=", Count)),
            size = 2.5, hjust = -0.1)+
  theme_bw()+
  ggforce::facet_col(facets = vars(Category), 
                     scales = "free_y", 
                     space = "free") +
  scale_fill_manual(values = c((viridis::plasma(4)[1]), (viridis::plasma(4)[4])))+
  xlab("-Log10(Adjusted P-value)")+
  ylab(NULL)+
  labs(fill = "Gene Ontology\n(GO) Categories")+
  theme(strip.text.x = element_text(angle = 0, hjust = 0, 
                                    size = 10, face = "bold"),
        strip.background = element_blank())+
  coord_cartesian(xlim=c(0,4)) +
  theme(plot.margin = margin(1.05,1.05,1.05,1.05, "cm"))

## NEGATIVE ----

## enrichment analysis using gene names
gost_results <- gost(query = model_genes_negative, 
                     organism = "drerio", 
                     multi_query = FALSE, 
                     evcodes = TRUE)

# modify the g:Profiler data frame
gost_mod_neg = gost_results$result[,c("query", "source", "term_id",
                                  "term_name", "p_value", "query_size", 
                                  "intersection_size", "term_size", 
                                  "effective_domain_size", "intersection")]

names(gost_mod_neg) = c("Cluster", "Category", "ID", "Description", "p.adjust",
                    "query_size", "Count", "term_size", "effective_domain_size",
                    "geneID")
gost_mod_neg$geneID = gsub(",", "/", gost_mod_neg$geneID)

## define as enrichResult object
gost_mod_neg_enrich  = new("enrichResult", result = gost_mod_neg)

## create plot data
plot_data <- gost_mod_neg_enrich@result %>%
  rownames_to_column(., var = "go_term") %>%
  mutate(neg_log10_padjust = -log10(p.adjust), 
         Category = gsub("GO:BP", "Biological process", Category) %>%
           gsub("GO:CC", "Cellular component", .) %>%
           gsub("GO:MF", "Molecular function", .)) %>%
  ## arrange by pvalue
  arrange(., desc(neg_log10_padjust)) %>%
  ## remove MIRNA
  filter(Category != "MIRNA")

## make category and ordered factor
plot_data$Description = factor(plot_data$Description, levels = plot_data$Description)

## plot as per ben's suggestion
negative_plot <- ggplot(plot_data, aes(x=neg_log10_padjust, 
                                       fill = Category,
                                       y=Description))+
  geom_bar(stat="identity")+
  geom_text(aes(label = paste0("n=", Count)),
            size = 2.5, hjust = -0.1)+
  theme_bw()+
  ggforce::facet_col(facets = vars(Category), 
                     scales = "free_y", 
                     space = "free") +
  scale_fill_manual(values = c((viridis::plasma(4)[1]), 
                               (viridis::plasma(4)[2]), 
                               (viridis::plasma(4)[3])))+
  xlab("-Log10(Adjusted P-value)")+
  ylab(NULL)+
  labs(fill = "Gene Ontology\n(GO) Categories")+
  theme(strip.text.x = element_text(angle = 0, hjust = 0,
                                    size = 10, face = "bold"),
        strip.background = element_blank())+
  coord_cartesian(xlim=c(0,max(plot_data$neg_log10_padjust*1.1))) +
  theme(plot.margin = margin(0.3,1.05,0.3,1.05, "cm"))

## combine plots
ggarrange(positive_plot+
            theme(plot.margin=margin(0.3, 1.05, 0.3, 4.2, "cm")),
          negative_plot,
          legend="none",
          nrow=2,
          heights=c(0.275, 0.725),
          labels=c("A", "B"))

## save them
ggsave("FiguresTables/Figure [go_terms_comb].tiff", 
       width = 21, height = 29.7/2, units="cm")

## save them
ggsave("FiguresTables/Figure [go_terms_comb].pdf", 
       width = 21, height = 29.7/2, units="cm")

## create supp table
gost_mod_pos$Weighting <- "positive"
gost_mod_neg$Weighting <- "negative"

supp <- full_join(gost_mod_pos, gost_mod_neg) %>%
  relocate(Weighting)

## export it
write.csv(supp, "FiguresTables/Table S[enrich].csv", row.names=FALSE)

## END SCRIPT