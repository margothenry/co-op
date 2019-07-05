library(ggplot2)
library(readr)
library(Hmisc)
library(RColorBrewer)
library(janitor)
library(here)
##############################

### Import analyses:
master_analyses <- read_csv(here("data_out", "master_analyses.csv"))
master_analyses$strain_info <- master_analyses$strain_info %>%
  replace(is.na(.), "0")
for (i in 1:nrow(master_analyses)) {
  if (master_analyses[i, "strain_info"] != "0"){
    master_analyses[i, "paper"] <- paste(master_analyses[i, "paper"], " ", "(", master_analyses[i, "strain_info"], ")", sep = "")
  }
}
# c-hyper vs generation (multiple generations):
### Most datasets use generations to notate timepoints. As of 190623, only 1 dataset uses "days", and 1 dataset uses "flasks". 
### Will need to clump all timepoints into more general ones (i.e. early, intermediate, late, etc.).
multiple_wide_generation_analysis <- master_analyses %>%
  subset(func == "multiple_wide") %>%
  subset(generation != "NA")

## A dot & line chart:
ggplot(multiple_wide_generation_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + geom_line(size = 0.75) + 
  scale_x_log10() + ggtitle("c-hyper vs generation") + theme(plot.title = element_text(hjust = 0.5, size = 11))


## Subset the generation-notating datasets by species:
### E. coli (both K-12 & ATCC):
multiple_wide_generation_Ecoli_analysis <- multiple_wide_generation_analysis %>%
  subset(species == "Ecoli_K12" | species == "Ecoli_ATCC")
ggplot(multiple_wide_generation_Ecoli_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + scale_x_log10() + 
  geom_line(size = 0.75) + ggtitle("c-hyper vs generation, E. coli") + theme(plot.title = element_text(hjust = 0.5, size = 11))
### Observation: Tendency towards a downtrend (exception: Sherlock2019, Wielgoss2016 (mucoid) (up); Kacar2017 (level))

### Sac:
multiple_wide_generation_Sac_analysis <- multiple_wide_generation_analysis %>%
  subset(species == "Sac")
ggplot(multiple_wide_generation_Sac_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + scale_x_log10() + 
  geom_line(size = 0.75) + ggtitle("c-hyper vs generation, S. cerevisiae") + theme(plot.title = element_text(hjust = 0.5, size = 11))
### Observation: There is a rising trend here but can't say much from just 2 datasets.


# Endpoints - c-hyper (?) - single gen, multi gen - superimpose on timeseries.