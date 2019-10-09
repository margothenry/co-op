library(ggplot2)
library(readr)
library(Hmisc)
library(RColorBrewer)
library(janitor)
library(here)
library (scales)
##############################
### (TL): Below are the codes I use to generate the visualizations in "data_out/images".
##############################
# Import analyses:
master_analyses <- read_csv(here("data_out", "master_analyses.csv"))
# Inspect the number of papers analyzed:
(paper_count <- length(unique(master_analyses$paper)))

# Inspect the number of papers, by species:
## E. coli:
Ecoli_entries_index <- grep("Ecoli", master_analyses$species)
Ecoli_entries <- master_analyses[Ecoli_entries_index,]
(Ecoli_paper_count <- length(unique(Ecoli_entries$paper)))

## S. cerevisiae:
Sac_entries_index <- grep("Sac", master_analyses$species)
Sac_entries <- master_analyses[Sac_entries_index,]
(Sac_paper_count <- length(unique(Sac_entries$paper)))

## P. aeruginosa:
aeruginosa_entries_index <- grep("aeruginosa", master_analyses$species)
aeruginosa_entries <- master_analyses[aeruginosa_entries_index,]
(aeruginosa_paper_count <- length(unique(aeruginosa_entries$paper)))

# Replace N/A values for "strain_info" with "0".
master_analyses$strain_info <- master_analyses$strain_info %>%
  replace(is.na(.), "0")

# If there are multiple strain backgrounds in 1 paper, include the value for "strain_info" in the "paper" column (for now).
for (i in 1:nrow(master_analyses)) {
  if (master_analyses[i, "strain_info"] != "0"){
    master_analyses[i, "paper"] <- paste(master_analyses[i, "paper"], " ", "(", master_analyses[i, "strain_info"], ")", sep = "")
  }
}

### Inspect the number of datasets (some papers have > 1 dataset):
(dataset_count <- length(unique(master_analyses$paper)))

# Inspect the number of datasets, by species:
## E. coli:
Ecoli_entries_strain_info_added <- master_analyses[Ecoli_entries_index,]
(Ecoli_dataset_count <- length(unique(Ecoli_entries_strain_info_added$paper)))

## S. cerevisiae:
Sac_entries_strain_info_added <- master_analyses[Sac_entries_index,]
(Sac_dataset_count <- length(unique(Sac_entries_strain_info_added$paper)))

## P. aeruginosa:
aeruginosa_entries_strain_info_added <- master_analyses[aeruginosa_entries_index,]
(aeruginosa_dataset_count <- length(unique(aeruginosa_entries_strain_info_added$paper)))

# c-hyper vs generation:
### Most datasets use generations to notate timepoints. However, there are still some that uses days or flasks. 
### Will need to clump all timepoints into more general ones (i.e. early/intermediate/late or convert everything to generations).
generation_analysis <- master_analyses %>%
  subset(generation != "NA")
## A dot & line chart, color coded by paper. Remove whatever elements you don't need:
ggplot(generation_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + labs(color = "Datasets") + geom_line(size = 0.75) + 
  scale_x_log10() + ggtitle("c-hyper vs generation") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), legend.position = "bottom") 
### A version optimized for the poster competition i.e. graph legend & title removed:
ggplot(generation_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + geom_line(size = 0.75) + 
  scale_x_log10() + theme(legend.position = "none")

## Split the chart by species:
### E. coli:
generation_analysis_Ecoli_index <- grep("E. coli", generation_analysis$species)
generation_analysis_Ecoli <- generation_analysis[generation_analysis_Ecoli_index,]
ggplot(generation_analysis_Ecoli, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + geom_line(size = 0.75) + 
  scale_x_log10() + theme(legend.position = "none")
### S. cerevisiae:
generation_analysis_Sac_index <- grep("S. cerevisiae", generation_analysis$species)
generation_analysis_Sac <- generation_analysis[generation_analysis_Sac_index,]
ggplot(generation_analysis_Sac, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + geom_line(size = 0.75) + 
  scale_x_log10() + theme(legend.position = "none")
### P. aeruginosa:
generation_analysis_P_aeruginosa_index <- grep("P. aeruginosa", generation_analysis$species)
generation_analysis_P_aeruginosa <- generation_analysis[generation_analysis_P_aeruginosa_index,]
ggplot(generation_analysis_P_aeruginosa, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + geom_line(size = 0.75) + 
  scale_x_log10() + theme(legend.position = "none")

# ## IN DEVELOPMENT - If you just need to color-code by species:
# ### First, we need to remove the strain names and just keep the general species name. 
# ### Also, use a neater-looking form of species names for the graph legend.
# #### E. coli:
# generation_analysis$species <- generation_analysis$species %>%
#   replace(grep("coli", generation_analysis$species), "E. coli")
# #### S. cerevisiae:
# generation_analysis$species <- generation_analysis$species %>%
#   replace(grep("Sac", generation_analysis$species), "S. cerevisiae")
# #### P. aeruginosa:
# generation_analysis$species <- generation_analysis$species %>%
#   replace(grep("aeruginosa", generation_analysis$species), "P. aeruginosa")
# ### Designate colors to each species:
# #### Since we designate color by species, why not set the stage for color-coding with the species?
# generation_analysis$graph_color <- generation_analysis$species
# #### Then, replace each species' name with the color you want:
# generation_analysis$graph_color <- generation_analysis$graph_color %>%
#   replace(grep("E. coli", generation_analysis$graph_color), "green")
# generation_analysis$graph_color <- generation_analysis$graph_color %>%
#   replace(grep("S. cerevisiae", generation_analysis$graph_color), "red")
# generation_analysis$graph_color <- generation_analysis$graph_color %>%
#   replace(grep("P. aeruginosa", generation_analysis$graph_color), "blue")
# ### Finally, the plot:
# ggplot(generation_analysis, aes(x = generation, y = c_hyper, color = graph_color)) + geom_point() + geom_line(size = 0.75) + 
#   labs(color = "Species") + scale_x_log10() + ggtitle("c-hyper vs generation") + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom")

# c-hyper vs generation (multiple generations):
multiple_wide_generation_analysis <- generation_analysis %>%
  subset(func == "multiple_wide")
ggplot(multiple_wide_generation_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + labs(color = "Datasets") + geom_line(size = 0.75) + 
  scale_x_log10() + ggtitle("c-hyper vs generation, experiments with multiple generations") + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom")
## Get the correlation between generation and c-hyper:
cor(multiple_wide_generation_analysis$c_hyper, multiple_wide_generation_analysis$generation)

# c-hyper vs generation (generation-notated end points only):
end_point_generation_analysis <- c()
end_point_generation_analysis <- generation_analysis %>%
  subset(func != "multiple_wide")
### Get only the end points of experiments with data for multiple generations (from multiple_wide_generation_analysis):
### [The selection of rows in this step is still manual & can be made automatic].
end_point_entries_multiple_wide <- c(6, 8, 12, 17, 28, 36, 46, 57, 59)
end_point_generation_analysis <- rbind(end_point_generation_analysis, multiple_wide_generation_analysis[end_point_entries_multiple_wide,])
### The plot is currently optimized for the poster competition (bigger axis.text, bigger legend.title, etc.)
ggplot(end_point_generation_analysis, aes(x = generation, y = c_hyper, color = species)) + xlab("Generations") + ylab("c-hyper") + 
  geom_jitter(size = 4, width = 0.2) + labs(color = "Species") + scale_x_log10() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18, face = "bold"), 
        legend.title = element_text(size = 18, face = "bold"), legend.position = "right", legend.text = element_text(size = 14, face = "italic")) 
# Include this if need a title: "+ ggtitle("c-hyper vs generation, generation-notating end points only") + theme(plot.title = element_text(hjust = 0.5, size = 11))" 
# Include this if need a linear regression line: "+ geom_smooth(method=lm, se = FALSE)"
# Include this if need a log curve: "+ geom_smooth(method = lm, formula = y ~ log(x), se = FALSE)"

## c-hyper vs generation (generation-notated end points only), by species:
### E. coli (also omitting Tenaillon2016 for the same reason as mentioned above):
end_point_generation_Ecoli_entries <- grep("E.", end_point_generation_analysis$species)
end_point_generation_Ecoli_analysis <- end_point_generation_analysis[end_point_generation_Ecoli_entries,]
ggplot(end_point_generation_Ecoli_analysis, aes(x = generation, y = c_hyper)) + geom_point() +  geom_smooth(method=lm, se=FALSE) + labs(color = "Datasets") + 
  ggtitle("c-hyper vs generation, generation-notating end points only, E. coli") + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom") 
### Sac:
end_point_generation_Sac_entries <- grep("Sac", end_point_generation_analysis$species)
end_point_generation_Sac_analysis <- end_point_generation_analysis[end_point_generation_Sac_entries,]
ggplot(end_point_generation_Sac_analysis, aes(x = generation, y = c_hyper)) + geom_point() +  geom_smooth(method=lm, se=FALSE) + labs(color = "Datasets") + 
  ggtitle("c-hyper vs generation, generation-notating end points only, S. cerevisiae") + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom") 
### P. aeruginosa:
end_point_generation_aeruginosa_entries <- grep("aeruginosa", end_point_generation_analysis$species)
end_point_generation_aeruginosa_analysis <- end_point_generation_analysis[end_point_generation_aeruginosa_entries,]
ggplot(end_point_generation_aeruginosa_analysis, aes(x = generation, y = c_hyper)) + geom_point() +  geom_smooth(method=lm, se=FALSE) + labs(color = "Datasets") + 
  ggtitle("c-hyper vs generation, generation-notating end points only, P. aeruginosa") + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom") 

