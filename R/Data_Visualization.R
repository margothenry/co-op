library(ggplot2)
library(readr)
library(Hmisc)
library(RColorBrewer)
library(janitor)
library(here)
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

master_analyses$strain_info <- master_analyses$strain_info %>%
  replace(is.na(.), "0")
### If there are multiple strain backgrounds in 1 paper, include the value for "strain_info" in the "paper" column (for now).
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

# c-hyper vs generation (all datasets):
## 190729 - FIX GRIFFITH2019 & JERISON2017 - DATASETS DIFFER IN SELECTIVE PRESSURE.
### Most datasets use generations to notate timepoints. However, there are still some that uses days or flasks. 
### Will need to clump all timepoints into more general ones (i.e. early/intermediate/late or convert everything to generations).
generation_analysis <- master_analyses %>%
  subset(generation != "NA")
## A dot & line chart:
ggplot(generation_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + labs(color = "Datasets") + geom_line(size = 0.75) + 
  scale_x_log10() + ggtitle("c-hyper vs generation") + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom") 

# c-hyper vs generation (multiple generations):
multiple_wide_generation_analysis <- generation_analysis %>%
  subset(func == "multiple_wide")
ggplot(multiple_wide_generation_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + labs(color = "Datasets") + geom_line(size = 0.75) + 
  scale_x_log10() + ggtitle("c-hyper vs generation, experiments with multiple generations") + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom")
cor(multiple_wide_generation_analysis$c_hyper, multiple_wide_generation_analysis$generation)

# c-hyper vs generation (generation-notated end points only):
end_point_generation_analysis <- generation_analysis %>%
  subset(func != "multiple_wide")
### Get only the end points of experiments with data for multiple generations:
### TL deliberately omitted the endpoint of Tenaillon2016 - way beyond the range of the other experiments, shifting the regression line too much.
### [The selection of rows in this step is still manual & can be made automatic].
end_point_entries_multiple_wide <- c(6, 8, 12, 17, 28, 36, 46, 59)
end_point_generation_analysis <- rbind(end_point_generation_analysis, multiple_wide_generation_analysis[end_point_entries_multiple_wide,])
ggplot(end_point_generation_analysis, aes(x = generation, y = c_hyper)) + geom_point() +  geom_smooth(method=lm, se=FALSE) + labs(color = "Datasets") + 
  ggtitle("c-hyper vs generation, generation-notating end points only") + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom") 

## End-point analysis, by species:
### E. coli (also omitting Tenaillon2016 for the same reason as mentioned above):
end_point_generation_Ecoli_entries <- grep("Ecoli", end_point_generation_analysis$species)
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

