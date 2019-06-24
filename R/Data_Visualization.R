library(ggplot2)
library(readr)
library(Hmisc)
library(RColorBrewer)
##############################

# c-hyper vs generation (multiple generations):
Tenaillon2016_analysis <- read_csv(here("data_out", "analyses", "Tenaillon2016_Analysis.csv"))
Lang2013_analysis <- read_csv(here("data_out", "analyses", "Lang2013_Analysis.csv"))
Sherlock2013_analysis <- read_csv(here("data_out", "analyses", "Sherlock2013_Analysis.csv"))
Sherlock2019_analysis <- read_csv(here("data_out", "analyses", "Sherlock2019_Analysis.csv"))
Tonoyan2019_analysis <- read_csv(here("data_out", "analyses", "Tonoyan2019_Analysis.csv"))
Wielgoss2016_mucoid_analysis <- read_csv(here("data_out", "analyses", "Wielgoss2016_mucoid_Analysis.csv"))
Wielgoss2016_nonMucoid_analysis <- read_csv(here("data_out", "analyses", "Wielgoss2016_nonMucoid_Analysis.csv"))
Sandberg2016_analysis <- read_csv(here("data_out", "analyses", "Sandberg2016_Analysis.csv"))
Kacar2017_analysis <- read_csv(here("data_out", "analyses", "Kacar2017_Analysis.csv"))
Morgenthaler2019_analysis <- read_csv(here("data_out", "analyses", "Morgenthaler2019_Analysis.csv"))

multiple_wide_analysis <- rbind(Tenaillon2016_analysis, Lang2013_analysis, Sherlock2013_analysis, Sherlock2019_analysis, Tonoyan2019_analysis,
                                Wielgoss2016_mucoid_analysis, Wielgoss2016_nonMucoid_analysis, Sandberg2016_analysis, Kacar2017_analysis, Morgenthaler2019_analysis)
### Most datasets use generations to notate timepoints. As of 190623, only 1 dataset uses "days", and 1 dataset uses "flasks". 
### Will need to clump all timepoints into more general ones (i.e. early, intermediate, late, etc.).
multiple_wide_generation_analysis <- multiple_wide_analysis %>%
  filter(generation != "NA")
myColors_generation <- brewer.pal(length(unique(multiple_wide_generation_analysis$paper)), "Set1")
names(myColors_generation) <- levels(multiple_wide_generation_analysis$paper)
colScale_generation <- scale_colour_manual(name = "paper", values = myColors_generation)

ggplot(multiple_wide_generation_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + 
  geom_line(size = 1.05) + colScale_generation + scale_x_log10() + ggtitle("c-hyper vs generation") + theme(plot.title = element_text(hjust = 0.5, size = 11))


## Subset the generation-notating datasets by species:
### E. coli (both K-12 & ATCC):
multiple_wide_generation_Ecoli_analysis <- filter(multiple_wide_generation_analysis, species != "Sac")
myColors_generation_Ecoli <- brewer.pal(length(unique(multiple_wide_generation_Ecoli_analysis$paper)), "Set1")
names(myColors_generation_Ecoli) <- levels(multiple_wide_generation_Ecoli_analysis$paper)
colScale_generation_Ecoli <- scale_colour_manual(name = "paper", values = myColors_generation_Ecoli)

multiple_wide_generation_Ecoli_analysis_graph <- ggplot(multiple_wide_generation_Ecoli_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + 
  geom_line(size = 1.05) + colScale_generation_Ecoli + scale_x_log10() + ggtitle("c-hyper vs generation, E. coli") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11))
### Observation: Tendency towards a downtrend (exception: Sherlock2019 (up), Kacar2017 (level))

### Sac:
multiple_wide_generation_Sac_analysis <- filter(multiple_wide_generation_analysis, species == "Sac")
myColors_generation_Sac <- brewer.pal(length(unique(multiple_wide_generation_Sac_analysis$paper)), "Set1")
names(myColors_generation_Sac) <- levels(multiple_wide_generation_Sac_analysis$paper)
colScale_generation_Sac <- scale_colour_manual(name = "paper", values = myColors_generation_Sac)

multiple_wide_generation_Sac_analysis_graph <- ggplot(multiple_wide_generation_Sac_analysis, aes(x = generation, y = c_hyper, color = paper)) + geom_point() + 
  geom_line(size = 1.05) + colScale_generation_Sac + scale_x_log10() + ggtitle("c-hyper vs generation, S. cerevisiae") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11))
### Observation: There is a rising trend here but can't say much from just 2 datasets.


# Endpoints - c-hyper (?) - single gen, multi gen:

# Endpoints - estimate (?)- single gen vs multi gen: