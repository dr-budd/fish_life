## script to create tree figure 

## set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## SET UP ----

## load libraries
# remotes::install_github("YuLab-SMU/ggtree")
# remotes::install_github("YuLab-SMU/ggtreeExtra")
library(ape) ## import tree
library(ggtree) ## plot tree
library(ggnewscale) ## enables second heatmap
library(viridis) ## colour blind friendly
library(ggtreeExtra) ## nice plots
library(tidyverse) ## helpful for data handling and visualization

## import summary metadata (don't use 09 because the order is not the same)
metadata_sum <- read.csv("dataFiles/06.00_metadata_sum.csv")

## import and format CpG density data
cpg_oe_data <- read.csv("dataFiles/06.00_cpg_oe_data.csv") %>%
  column_to_rownames(., var="X")

## quick sanity check
all(rownames(cpg_oe_data) == unique(metadata_sum$organism_name))

## import predictive promoters
predictive_promoters <- read.csv("dataFiles/10.00_predictive_promoters_avg.csv")

## FORMAT DATA AND CREATE TREE ----

## create list
pp_list <- predictive_promoters %>% pull(promoter_id)

## subset cg density data for list
cpg_oe_less <- cpg_oe_data[,pp_list]

## format long
cpg_oe_df <- cpg_oe_less %>%
  rownames_to_column(., var = "organism_name") %>%
  pivot_longer(., cols = -organism_name, 
               names_to = "promoter_id", 
               values_to = "cpg_oe")

## get timetree data
chordate_tree <- read.tree("dataFiles/chordata_species.nwk")

## edit tip labels
chordate_tree$tip.label <- gsub("_", " ", chordate_tree$tip.label)

## keep only intersecting tips
keep_tips <- intersect(chordate_tree$tip.label, metadata_sum$organism_name)

## subset tree for genome and lifespan info species
my_tree <- keep.tip(chordate_tree, keep_tips)

## remove chordate tree (it's big)
rm(chordate_tree)

## print tree
ggtree(my_tree)+
  geom_tiplab(size = 1)

## print rooted tree
ggtree(root(my_tree, 
            outgroup = "Danio rerio"))+
  geom_tiplab(size = 1)

## add node labels to check clade positions
ggtree(root(my_tree,
            outgroup = "Danio rerio"),
       size = 0.1)+ ## reduce branch thickness
  geom_tiplab(size = 1)+
  geom_text(aes(label=node), size=2)

## flip tree so that it is ordered by distance
## from zebrafish (because it will be presented
## as a chronogram...)
flip(ggtree(root(my_tree, 
                 outgroup = "Danio rerio")), 
     442, 444)

## fish mentioned in publication (to create tip labels)
pub_fish <- c("Danio rerio", 
              "Nothobranchius furzeri",
              "Sebastes aleutianus", 
              "Hoplostethus atlanticus", 
              "Danio kyathit",
              "Lethenteron camtschaticum",
              "Paralichthys olivaceus", 
              "Noturus placidus")

## make a semi-circular tree (with flipped node)
c <- flip(ggtree(root(my_tree, 
                 outgroup = "Danio rerio"), ## (root to zebra)
            size = 0.1,
            ## make it a semicircle
            layout = "fan",
            open.angle=180,
            ## draw cladogram
            branch.length="none"), 
          442, 444)

## print tree
# c

## format metdata to add to tree (for the heatmap and bar plot)
metadata_sum_matching <- metadata_sum %>%
  filter(organism_name %in% my_tree$tip.label) %>%
  relocate(organism_name) %>% ## IMPORTANT!
  group_by(order) %>%
  add_count(.) %>%
  # mutate(order_more = ifelse(n > 10, 
  #                            order, 
  #                            "other")) %>%
  # mutate(genus_s = ifelse(genus == "Sebastes", 
  #                         genus, 
  #                         "Other")) %>%
  # ungroup(.) %>%
  # group_by(genus) %>%
  # add_count(., name = "genus_count") %>%
  # mutate(genus_more = ifelse(genus_count > 10, 
  #                            genus, 
  #                            "Other")) %>%
  mutate(pub_fish = ifelse(organism_name %in% pub_fish, 
                           organism_name, 
                           NA))

## add formatted metadata to tree
c2 <- c %<+% metadata_sum_matching

## print tree
# c2

## for sebastes

# ## get viridis colours you want (drop the first colour)
# viridis_colours <- viridis(option = "D", n = 3)[2:3]
# 
# c3 <- c2+geom_tippoint(mapping=aes(x=x+0.75, colour=genus_s),
#                        size=0.75,
#                        alpha = 0.8,
#                        na.rm=TRUE,
#                        stroke=0)+
#   scale_colour_manual("Genus",
#                       values = viridis_colours,
#                       guide=guide_legend(keywidth=0.3,
#                                          keyheight=0.3,
#                                          ncol=1,
#                                          size = 1,
#                                          override.aes=list(size=2,alpha=1),
#                                          order=1))
# 
# c3

## for fish mentioned in the pub

## create vector of viridis colours (drop the first colour)
viridis_colours <- viridis(option = "turbo", n = 9)[2:9]

## create tree with coloured points for pub_fish
## (to help with labelling)
c3_temp <- c2+geom_tippoint(mapping=aes(x=x+0.75, 
                                        colour=pub_fish),
                       size=0.75,
                       alpha = 0.8,
                       stroke=0)+
  scale_colour_manual("Fish mentioned in-text",
                      values = viridis_colours,
                      na.value = "white",
                      guide=guide_legend(keywidth=0.3,
                                         keyheight=0.3,
                                         ncol=1,
                                         size=1,
                                         override.aes=list(size=2,alpha=1),
                                         order=1))

## print tree
c3_temp

## create tree with circles for pub_fish
## (for the actual figure)
c3 <- c2+geom_tippoint(mapping=aes(x=x+0.75, size=pub_fish),
                       # size=0.75,
                       alpha = 0.8,
                       shape=21,
                       stroke=0.1)+
  scale_size_manual("Fish mentioned in-text",
                    # values = c("black, black"),
                    # na.value = "white",
                    values=c(rep(0.75, 9), 0),
                    na.value = 0,
                    guide="none")

## print tree
# c3

## add CG density and lifespan metadata
## (heatmap and barplot)
c4 <- c3+
  ## add promoter CG density
  geom_fruit(data=cpg_oe_df, 
             geom=geom_tile,
             mapping=aes(y=organism_name, 
                         x=promoter_id, 
                         fill=cpg_oe),
             offset = 0.025, 
             pwidth = 0.4)+
  scale_fill_viridis_c(option = "C", 
                       breaks=c(0, 0.5, 1, 1.5, 2),
                       labels=c("0", "0.5", "1", "1.5", "> 2"),
                       na.value="transparent",
                       limits=c(0, 2),
                       "CpG O/E")+
  ## add lifespan bars
  geom_fruit(geom=geom_bar,
             mapping=aes(y=organism_name, 
                         x=mean_lifespan),
             fill="#0D0887FF",
             pwidth=0.38, 
             orientation="y", 
             stat="identity")+
  ## adjust the legend and margin
  theme(legend.position=c(0.875, 0.7),
        legend.key.size = unit(0.2, 'cm'),
        legend.title=element_text(size=5),
        legend.text=element_text(size=4), 
        plot.margin = unit(c(0,0,0,0), "cm"))

# ## print figure (inefficient)
# Sys.time()
# c4
# Sys.time()

## save it (takes a while)
print("save time: ")
Sys.time()
ggsave("figuresTables/Figure [tree].pdf", c4, width = 21, height = 29.7/2, units = "cm", limitsize=FALSE)
Sys.time()

## crop it
print("crop time: ")
Sys.time()
system2(command = "pdfcrop", args = c("Figure [tree].pdf","Figure [tree].pdf","--margins","10"))
knitr::plot_crop("figuresTables/Figure [tree].pdf")
Sys.time()

## to fix colours:

# ## plot cpgoe distribution
# plot(cpg_oe_df$cpg_oe)
# 
# ## so most things are 2.5 or less? 
# 
# plot(cpg_oe_df %>%
#        filter(cpg_oe > 2) %>%
#        pull(cpg_oe))

## END SCRIPT
