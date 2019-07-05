install.packages("Hmisc")
install.packages("plyr")
install.packages("naniar")
library(here)
library(dgconstraint)
library(tidyverse)
library(R.utils)
library(janitor)
library(tidyverse)
library(readr)
library(devtools)
library(Hmisc)
library(naniar)
source("R/dgconstraint/functions/multiple_wide.R")
source("R/dgconstraint/functions/multiple_long.R")
source("R/dgconstraint/functions/single_wide.R")
source("R/dgconstraint/functions/single_long.R")

############################
# (TL): Note: All files imported for the analysis must have "_usable.csv" at the end, for uniformity.
############################

############################
# Multiple wide:
############################
### (TL): A template for func calling (use EITHER generations, days, or flasks; use strain_info if paper has datasets from multiple founder strains):
### (TL): multiple_wide(paper = "", dataset_name = "", environment = "", selective_pressure = "", species = "", ploidy = "", (strain_info = ""), generations/days/flasks = c(""))

# (Working) Tenaillon2016:
### This paper originally had two clones sequenced for each timepoint. Collapse the two columns and save the resulting dataset.
Tenaillon2016 <- read_csv(here("data_in", "original & usable", "Tenaillon2016", "Tenaillon2016_usable.csv"))
Tenaillon2016 <- clean_names(Tenaillon2016, case = "snake")
colnames(Tenaillon2016) <- tolower(colnames(Tenaillon2016))
Tenaillon2016 <- Tenaillon2016 %>% 
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic")
### (TL): Why there's a second pipe here: If there are any cols removed from filter(), those changes need to be saved in the dataset before transmute(), which grabs cols using "$".
### (TL): In this particular case though, there are cols that need collapsing, so no direct col-grabbing w/ "$". You'll see examples of that eventually.
### (TL): If there are cols that need collapsing or renaming, use transmute(). If not, use select() - much less of a hassle.
Tenaillon2016 <- Tenaillon2016 %>% 
  transmute(gene, population, "500" = `x500_i1_r1`+`x500_i2_r1`, "1000" =`x1000_i1_r1`+`x1000_i2_r1`, "1500" =`x1500_i1_r1`+`x1500_i2_r1`,  "2000"= `x2000_i1_r1`+`x2000_i2_r1`, "5000"=`x5000_i1_r1`+`x5000_i2_r1`, "10000"= `x10000_i1_r1`+`x10000_i2_r1`, "15000"=`x15000_i1_r1`+`x15000_i2_r1`,"20000"=`x20000_i1_r1`+`x20000_i2_r1`,"30000"= `x30000_i1_r1`+`x30000_i2_r1`,"40000"=`x40000_i1_r1`+`x40000_i2_r1`,"50000"=`x50000_i1_r1`+`x50000_i2_r1`)
### (TL): Create a new variable to store the index of rows that contain intergenic mutations not notated by breseq.
### (TL): If there are any rows of this kind, remove them from the dataset. The if() was used because the dataset will be empty if there are no badly-notated rows (dataset(-0,) is empty!)
### (TL): grep() is used for all hitherto known patterns indicating intergenicness in the colum "gene". Future updates could include more patterns, depending on how creative people could get with their notations.
Tenaillon2016_out <- c(grep(",", Tenaillon2016$gene), grep("/", Tenaillon2016$gene))
if (length(Tenaillon2016_out) > 0) {
  Tenaillon2016 <- Tenaillon2016[-Tenaillon2016_out,]
}

write_csv(Tenaillon2016, here("data_in", "for_func", "Tenaillon2016.csv"))
### (TL): Make sure to call out ALL arguments - there are lots of parameter. so position-based value input isn't reliable.
multiple_wide(paper = "Tenaillon2016", dataset_name = "Tenaillon2016", environment = "Davis minimal medium", 
              generations = c("500", "1000", "1500", "2000", '5000', '10000', '15000', '20000', '30000', '40000','50000'), selective_pressure = "Davis minimal medium", 
              species = "Ecoli_K12", ploidy = "haploid")



# (Working) Lang2013:
Lang2013 <- read_csv(here("data_in", "original & usable", "Lang2013", "Lang2013_usable.csv"))
Lang2013 <- clean_names(Lang2013, case = "snake")
colnames(Lang2013) <- tolower(colnames(Lang2013))
### (TL): (Custom) If there are too many non-numeric timepoints to rename manually, gsub() for the naming pattern of a particular dataset, replacing the non-numbers with "" (no space).
names(Lang2013) <- gsub("p1_", "", names(Lang2013))
Lang2013 <- Lang2013 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic", gene != "0")
Lang2013 <- Lang2013 %>%
  select(gene, population, '0', '140', '240', '335', '415', '505', '585', '665', '745', '825', '910', '1000')
Lang2013_out <- c(grep(",", Lang2013$gene), grep("/", Lang2013$gene))
if (length(Lang2013_out) > 0) {
  Lang2013 <- Lang2013[-Lang2013_out,]
}

write_csv(Lang2013, here("data_in", "for_func", "Lang2013.csv"))
multiple_wide(paper = "Lang2013", dataset_name = "Lang2013", environment = "YPD", generations = c('0', '140', '240', '335', '415', '505', '585', '665', '745', '825', '910', '1000'), 
              selective_pressure = "YPD", species = "Sac", ploidy = "haploid")



# (Working) Sherlock2013:
Sherlock2013 <- read_csv(here("data_in", "original & usable", "Sherlock2013", "Sherlock2013_usable.csv"))
Sherlock2013 <- clean_names(Sherlock2013, case = "snake")
colnames(Sherlock2013) <- tolower(colnames(Sherlock2013))
### (TL): The renaming of non-numeric timepoint names is a bit TLcky here. Both "gene" & timepoint names have the letter "g" (both lowercase after tolower()).
### (TL): So, it's necessary to capitalize the "g" in "gene" so that we could gsub() the non-numeric timepoint names:
names(Sherlock2013) <- gsub("gene", "Gene", names(Sherlock2013))
names(Sherlock2013) <- gsub("g", "", names(Sherlock2013))
### (TL): Now, change the "Gene" column back to lowercase:
names(Sherlock2013) <- gsub("Gene", "gene", names(Sherlock2013))
Sherlock2013 <- Sherlock2013 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") 
Sherlock2013 <- Sherlock2013 %>%
  select(gene, population, "7", "70", "133", "196", "266", "322", "385", "448")
Sherlock2013_out <- c(grep(",", Sherlock2013$gene), grep("/", Sherlock2013$gene))
if (length(Sherlock2013_out) > 0) {
  Sherlock2013 <- Sherlock2013[-Sherlock2013_out,]
}

write_csv(Sherlock2013, here("data_in", "for_func", "Sherlock2013.csv"))
multiple_wide(paper = "Sherlock2013", dataset_name = "Sherlock2013", environment = "YPD", generations = c("7", "70", "133", "196", "266", "322", "385", "448"), 
              selective_pressure = "YPD", species = "Sac", ploidy = "haploid")



# (Working) Sherlock2019:
Sherlock2019 <- read_csv(here("data_in", "original & usable", "Sherlock2019", "Sherlock2019_usable.csv"))
Sherlock2019 <- clean_names(Sherlock2019, case = "snake")
colnames(Sherlock2019) <- tolower(colnames(Sherlock2019))
### (TL): If the input dataset already has numeric timepoint col names, clean_name() will add an "x" to each of them. Use gsub() to remove this.
names(Sherlock2019) <- gsub("x", "", names(Sherlock2019))
Sherlock2019 <- Sherlock2019 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic", gene != "0")
Sherlock2019 <- Sherlock2019 %>%
  select(gene, population, "0", "50", "100", "150", "200", "250", "300", "350", "400", "450", "500")
Sherlock2019_out <- c(grep(",", Sherlock2019$gene), grep("/", Sherlock2019$gene))
if (length(Sherlock2019_out) > 0) {   
  Sherlock2019 <- Sherlock2019[-Sherlock2019_out,] 
} 
write_csv(Sherlock2019, here("data_in", "for_func", "Sherlock2019.csv"))
multiple_wide(paper = "Sherlock2019", dataset_name = "Sherlock2019", environment = "Davis minimal medium", 
              generations = c('0', "50", "100", "150", "200", "250", "300", "350", "400", "450", "500"), selective_pressure = "Davis minimal medium", 
              species = "Ecoli_K12", ploidy = "haploid")



# (Working) Tonoyan2019:
Tonoyan2019 <- read_csv(here("data_in", "original & usable", "Tonoyan2019", "Tonoyan2019_usable.csv"))
Tonoyan2019 <- clean_names(Tonoyan2019, case = "snake")
colnames(Tonoyan2019) <- tolower(colnames(Tonoyan2019))
names(Tonoyan2019) <- gsub("x", "", names(Tonoyan2019))
Tonoyan2019 <- Tonoyan2019 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") 
Tonoyan2019 <- Tonoyan2019 %>%  
  select(gene, population, "1", "14", "20")
Tonoyan2019_out <- c(grep(",", Tonoyan2019$gene), grep("/", Tonoyan2019$gene))
if (length(Tonoyan2019_out) > 0) {
  Tonoyan2019 <- Tonoyan2019[-Tonoyan2019_out,] 
} 

write_csv(Tonoyan2019, here("data_in", "for_func", "Tonoyan2019.csv"))
multiple_wide(paper = "Tonoyan2019", dataset_name = "Tonoyan2019", environment = "LB medium", generations = c("1", "14", "20"), selective_pressure = "LB medium", 
              species = "Ecoli_ATCC", ploidy = "haploid")



# (Working) Wielgoss2016, subsetted by presence of mucoid:
## (Working) Wielgoss2016_mucoid:
Wielgoss2016_mucoid <- read_csv(here("data_in", "original & usable", "Wielgoss2016", "Wielgoss2016_mucoid_usable.csv"))
Wielgoss2016_mucoid <- clean_names(Wielgoss2016_mucoid, case = "snake")
colnames(Wielgoss2016_mucoid) <- tolower(colnames(Wielgoss2016_mucoid))
names(Wielgoss2016_mucoid) <- gsub("x", "", names(Wielgoss2016_mucoid))
Wielgoss2016_mucoid <- Wielgoss2016_mucoid %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic")
Wielgoss2016_mucoid <- Wielgoss2016_mucoid %>%
  select(gene, population, "0", "10")
Wielgoss2016_mucoid_out <- c(grep(",", Wielgoss2016_mucoid$gene), grep("/", Wielgoss2016_mucoid$gene))
if (length(Wielgoss2016_mucoid_out) > 0) {   
  Wielgoss2016_mucoid <- Wielgoss2016_mucoid[-Wielgoss2016_mucoid_out,] 
} 

write_csv(Wielgoss2016_mucoid, here("data_in", "for_func", "Wielgoss2016_mucoid.csv"))
multiple_wide(paper = "Wielgoss2016", dataset_name = "Wielgoss2016_mucoid", environment = "Lysogeny broth", generations = c("0", "10"), 
              selective_pressure = "Lysogeny broth", species = "Ecoli_K12", ploidy = "haploid", strain_info = "mucoid")


## (Working) Wielgoss2016_nonMucoid:
Wielgoss2016_nonMucoid <- read_csv(here("data_in", "original & usable", "Wielgoss2016", "Wielgoss2016_nonMucoid_usable.csv"))
Wielgoss2016_nonMucoid <- clean_names(Wielgoss2016_nonMucoid, case = "snake")
colnames(Wielgoss2016_nonMucoid) <- tolower(colnames(Wielgoss2016_nonMucoid))
names(Wielgoss2016_nonMucoid) <- gsub("x", "", names(Wielgoss2016_nonMucoid))
Wielgoss2016_nonMucoid <- Wielgoss2016_nonMucoid %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic")
Wielgoss2016_nonMucoid <- Wielgoss2016_nonMucoid %>%
  select(gene, population, "0", "10")
Wielgoss2016_nonMucoid_out <- c(grep(",", Wielgoss2016_nonMucoid$gene), grep("/", Wielgoss2016_nonMucoid$gene))
if (length(Wielgoss2016_nonMucoid_out) > 0) {   
  Wielgoss2016_nonMucoid <- Wielgoss2016_nonMucoid[-Wielgoss2016_nonMucoid_out,] 
} 

write_csv(Wielgoss2016_nonMucoid, here("data_in", "for_func", "Wielgoss2016_nonMucoid.csv"))
multiple_wide(paper = "Wielgoss2016", dataset_name = "Wielgoss2016_nonMucoid", environment = "Lysogeny broth", generations = c("0", "10"), 
              selective_pressure = "Lysogeny broth", species = "Ecoli_K12", ploidy = "haploid", strain_info = "non-mucoid")



# (Working) Sandberg2016:
Sandberg2016 <- read_csv(here("data_in", "original & usable", "Sandberg2016", "Sandberg2016_usable.csv"))
Sandberg2016 <- clean_names(Sandberg2016, case = "snake")
colnames(Sandberg2016) <- tolower(colnames(Sandberg2016))
Sandberg2016 <- Sandberg2016 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic")
Sandberg2016 <- Sandberg2016 %>%
  transmute(gene, population, "23" = Sandberg2016$"flask_23", "58" = Sandberg2016$"flask_58", "133" = Sandberg2016$"flask_133")
Sandberg2016_out <- c(grep(",", Sandberg2016$gene), grep("/", Sandberg2016$gene))
if (length(Sandberg2016_out) > 0) {   
  Sandberg2016 <- Sandberg2016[-Sandberg2016_out,] 
} 

write_csv(Sandberg2016, here("data_in", "for_func", "Sandberg2016.csv"))
### (TL): Make sure to call out the argument "flasks"!
multiple_wide(paper = "Sandberg2016", dataset_name = "Sandberg2016", environment = "Davis minimal medium", selective_pressure = "Davis minimal medium", 
              species = "Ecoli_K12", ploidy = "haploid", flasks = c("23", "58", "133"))



# (TL): Kacar2017:
Kacar2017 <- read_csv(here("data_in", "original & usable", "Kacar2017", "Kacar2017_usable.csv"))
Kacar2017 <- clean_names(Kacar2017, case = "snake")
colnames(Kacar2017) <- tolower(colnames(Kacar2017))
names(Kacar2017) <- gsub("x", "", names(Kacar2017))
Kacar2017 <- Kacar2017 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic", gene != "0")
Kacar2017 <- Kacar2017 %>%
  select(gene, population, "500", "1000", "1500", "2000")
Kacar2017_out <- c(grep(",", Kacar2017$gene), grep("/", Kacar2017$gene))
if (length(Kacar2017_out) > 0) {   
  Kacar2017 <- Kacar2017[-Kacar2017_out,] 
  } 

write_csv(Kacar2017, here("data_in", "for_func", "Kacar2017.csv"))
multiple_wide(paper = "Kacar2017", dataset_name = "Kacar2017", environment = "Minimal glucose medium", generations = c("500", "1000", "1500", "2000"), 
              selective_pressure = "Minimal glucose medium", species = "Ecoli_K12", ploidy = "haploid")



# (TL): Morgenthaler2019:
Morgenthaler2019 <- read_csv(here("data_in", "original & usable", "Morgenthaler2019", "Morgenthaler2019_usable.csv"))
Morgenthaler2019 <- clean_names(Morgenthaler2019, case = "snake")
colnames(Morgenthaler2019) <- tolower(colnames(Morgenthaler2019))
Morgenthaler2019 <- Morgenthaler2019 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic")
Morgenthaler2019 <- Morgenthaler2019 %>%
  transmute(gene, population, "42" = "day_42", "50" = "day_50")
Morgenthaler2019_out <- c(grep(",", Morgenthaler2019$gene), grep("/", Morgenthaler2019$gene))
if (length(Morgenthaler2019_out) > 0) {   
  Morgenthaler2019 <- Morgenthaler2019[-Morgenthaler2019_out,] 
} 

write_csv(Morgenthaler2019, here("data_in", "for_func", "Morgenthaler2019.csv"))
multiple_wide(paper = "Morgenthaler2019", dataset_name = "Morgenthaler2019", environment = "M9", days = c("42", "50"), selective_pressure = "M9", 
              species = "Ecoli_K12", ploidy = "haploid")



# (TL): Du2019:
Du2019 <- read_csv(here("data_in", "original & usable", "Du2019", "Du2019_usable.csv"))
Du2019 <- clean_names(Du2019, case = "snake")
colnames(Du2019) <- tolower(colnames(Du2019))
names(Du2019) <- gsub("x", "", names(Du2019))
Du2019 <- Du2019 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic", gene != "0")
### (TL): To assist with the spread of the "flask" column, I manually added a "value" column in Excel. All values for the "value" column are "1".
### (TL): Spread the "flask" column:
Du2019 <- spread(Du2019, flask, value)
### (TL): Gather the population columns ("aa1" to "aa6"):
Du2019 <- gather(Du2019, population, frequency, "aa1" : "aa6", factor_key=TRUE)
Du2019 <- Du2019 %>%
  replace(is.na(.), 0)
Du2019 <- Du2019 %>%
  select(gene, population, "intermediate", "late")
Du2019_out <- c(grep(",", Du2019$gene), grep("/", Du2019$gene))
if (length(Du2019_out) > 0) {   
  Du2019 <- Du2019[-Du2019_out,] 
} 

write_csv(Du2019, here("data_in", "for_func", "Du2019.csv"))
multiple_wide(paper = "Du2019", dataset_name = "Du2019", environment = "Minimal glucose medium", flasks = c("intermediate", "late"), 
              selective_pressure = "pH 5.5", species = "Ecoli_K12", ploidy = "haploid")



# (TL): Flynn2014, subset by founder:
## (TL): Flynn2014_biofilm:
Flynn2014_biofilm <- read_csv(here("data_in", "original & usable", "Flynn2014", "Flynn2014_biofilm_usable.csv"))
Flynn2014_biofilm <- clean_names(Flynn2014_biofilm, case = "snake")
colnames(Flynn2014_biofilm) <- tolower(colnames(Flynn2014_biofilm))
names(Flynn2014_biofilm) <- gsub("x", "", names(Flynn2014_biofilm))
Flynn2014_biofilm <- Flynn2014_biofilm %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic", gene != "0")
Flynn2014_biofilm <- Flynn2014_biofilm %>%
  select(gene, population, "102", "150", "264", "396", "450", "540")
Flynn2014_biofilm_out <- c(grep(",", Flynn2014_biofilm$gene), grep("/", Flynn2014_biofilm$gene))
if (length(Flynn2014_biofilm_out) > 0) {   
  Flynn2014_biofilm <- Flynn2014_biofilm[-Flynn2014_biofilm_out,] 
} 

write_csv(Flynn2014_biofilm, here("data_in", "for_func", "Flynn2014_biofilm.csv"))
multiple_wide(paper = "Flynn2014", dataset_name = "Flynn2014_biofilm", environment = "M63", generations = c("102", "150", "264", "396", "450", "540"), 
              selective_pressure = "Polystyrene beads", species = "P_aeruginosa_PA14", ploidy = "haploid", strain_info = "biofilm-evolved")


## (TL): Flynn2014_planktonic:
Flynn2014_planktonic <- read_csv(here("data_in", "original & usable", "Flynn2014", "Flynn2014_planktonic_usable.csv"))
Flynn2014_planktonic <- clean_names(Flynn2014_planktonic, case = "snake")
colnames(Flynn2014_planktonic) <- tolower(colnames(Flynn2014_planktonic))
names(Flynn2014_planktonic) <- gsub("x", "", names(Flynn2014_planktonic))
Flynn2014_planktonic <- Flynn2014_planktonic %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic", gene != "0")
Flynn2014_planktonic <- Flynn2014_planktonic %>%
  select(gene, population, "396", "540")
Flynn2014_planktonic_out <- c(grep(",", Flynn2014_planktonic$gene), grep("/", Flynn2014_planktonic$gene))
if (length(Flynn2014_planktonic_out) > 0) {   
  Flynn2014_planktonic <- Flynn2014_planktonic[-Flynn2014_planktonic_out,] 
} 

write_csv(Flynn2014_planktonic, here("data_in", "for_func", "Flynn2014_planktonic.csv"))
multiple_wide(paper = "Flynn2014", dataset_name = "Flynn2014_planktonic", environment = "M63", generations = c("396", "540"), 
              selective_pressure = "No polystyrene beads", species = "P_aeruginosa_PA14", ploidy = "haploid", strain_info = "planktonic-evolved")



# (TL): Keane2014:
Keane2014 <- read_csv(here("data_in", "original & usable", "Keane2014", "Keane2014_usable.csv"))
Keane2014 <- clean_names(Keane2014, case = "snake")
colnames(Keane2014) <- tolower(colnames(Keane2014))
names(Keane2014) <- gsub("x", "", names(Keane2014))
### (TL): Replace "0" values with "NA" - many funcs are built for NA's, so it's easier to mass-replace NA's than some other values:
Keane2014 <- Keane2014 %>%
  replace_with_na(replace = list("0" = 0, "440" = 0, "1100" = 0, "1540" = 0, "1980" = 0, "2200" = 0)) %>% 
  filter(details != "intergenic", details != "0")
### (TL): Change non-zero values in pop cols to "1", then change NA's back to 0 to feed in the func.
Keane2014[, 4:ncol(Keane2014)] <- Keane2014[, 4:ncol(Keane2014)] %>%
  replace(!is.na(.), 1) %>%
  replace(is.na(.), 0)
Keane2014 <- Keane2014 %>%
  select(gene, population, "0", "440", "1100", "1540", "1980", "2200")
Du2019_out <- c(grep(",", Du2019$gene), grep("/", Du2019$gene))
if (length(Du2019_out) > 0) {   
  Du2019 <- Du2019[-Du2019_out,] 
} 

write_csv(Du2019, here("data_in", "for_func", "Du2019.csv"))
multiple_wide(paper = "Du2019", dataset_name = "Du2019", environment = "Minimal glucose medium", flasks = c("intermediate", "late"), 
              selective_pressure = "pH 5.5", species = "Ecoli_K12", ploidy = "haploid")

#####################################
# Multiple long:
#####################################
### (TL): A template for func calling:
### (TL): multiple_long(paper = "", dataset_name = "", environment = "", generations/days/flasks = "", selective_pressure = c(""), species = "", ploidy = "", (strain_info = ""))

# (Working) Jerison2017:
Jerison2017 <- read_csv(here("data_in", "original & usable", "Jerison2017", "Jerison2017_usable.csv"))
Jerison2017 <- clean_names(Jerison2017, case = "snake")
colnames(Jerison2017) <- tolower(colnames(Jerison2017))
Jerison2017 <- Jerison2017 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic", gene != "0") %>%
  filter(frequency != "clone")
Jerison2017 <- Jerison2017 %>%
  select(gene, population, selective_pressure, frequency)
Jerison2017_out <- c(grep(",", Jerison2017$gene), grep("/", Jerison2017$gene))
if (length(Jerison2017_out) > 0) {   
  Jerison2017 <- Jerison2017[-Jerison2017_out,] 
} 

write_csv(Jerison2017, here("data_in", "for_func", "Jerison2017.csv"))
multiple_long(paper = "Jerison2017", dataset_name = "Jerison2017", environment = "SC", generations = "500", selective_pressure = c("HT", "OT"), 
              species = "Sac", ploidy = "haploid")



# (Working) Payen2016, subset by ploidy:
Payen2016 <- read_csv(here("data_in", "original & usable", "Payen2016", "Payen2016_usable.csv"))
Payen2016 <- clean_names(Payen2016, case = "snake")
colnames(Payen2016) <- tolower(colnames(Payen2016))
Payen2016 <- Payen2016 %>%
  replace(is.na(.), 0) %>% 
  ### (TL): Besides the usual "intergenic" screen, filter out the values whose frequencies are "clone" & the values w/o a gene.
  ### (TL): Also, this particular dataset is a special case - all the frequencies of the "sulfate" selective_pressure are "0", which would lead to an error if put through the func.
  ### (TL): To fix this, we need to also filter out the values whose frequencies are "0" [Not sure if this should be included in all analyses - theoretically speaking, this filter would make the func run faster].
  filter(details != "intergenic", frequency != "clone", frequency != "0", gene != "0")
  
Payen2016_out <- c(grep(",", Payen2016$gene), grep("/", Payen2016$gene))
if (length(Payen2016_out) > 0) {   
  Payen2016 <- Payen2016[-Payen2016_out,] 
}

## Payen2016_dip:
Payen2016_dip <- Payen2016 %>%
  ### (TL): Filter for diploid observations:
  filter(ploidy == "diploid") %>%
  select(gene, population, selective_pressure, frequency)

write_csv(Payen2016_dip, here("data_in", "for_func", "Payen2016_dip.csv"))
multiple_long(paper = "Payen2016", dataset_name = "Payen2016_dip", environment = "phosphate", generations = "20", selective_pressure = c("phosphate"), species = "Sac", ploidy = "diploid")


## Payen2016_hap:
Payen2016_hap <- Payen2016 %>%
  ### (TL): Filter for haploid observations:
  filter(ploidy == "haploid") %>%
  select(gene, population, selective_pressure, frequency)

write_csv(Payen2016_hap, here("data_in", "for_func", "Payen2016_hap.csv"))
multiple_long(paper = "Payen2016", dataset_name = "Payen2016_hap", environment = "LB", generations = "20", selective_pressure = c("phosphate"), species = "Sac", ploidy = "haploid")



# (TL): (Working) Lennen2015, subsetted by founder (founders have different plasmids):
## (TL): Lennen2015_pMA1:
Lennen2015_pMA1 <- read_csv(here("data_in", "original & usable", "Lennen2015", "Lennen2015_pMA1_usable.csv"))
Lennen2015_pMA1 <- clean_names(Lennen2015_pMA1, case = "snake")
colnames(Lennen2015_pMA1) <- tolower(colnames(Lennen2015_pMA1))
### (TL): Convert to multiple_long:
Lennen2015_pMA1 <- gather(Lennen2015_pMA1, population, frequency, "m9p_ma1_1" : "na_clp_ma1_19", factor_key=TRUE)
### (TL): Add a new col, "selective_pressure". The value "0" is just a placeholder:
Lennen2015_pMA1 <- cbind(Lennen2015_pMA1, selective_pressure = 0)
### (TL): Get the rows that have M9 in the population name & save their values in "selective_pressure" as "no NaCl". Same thing for NaCl:
Lennen2015_pMA1_no_NaCl_row_nums <- grep("m9", Lennen2015_pMA1$population)
Lennen2015_pMA1[Lennen2015_pMA1_no_NaCl_row_nums, "selective_pressure"] <- "no NaCl"
Lennen2015_pMA1_NaCl_row_nums <- grep("na_cl", Lennen2015_pMA1$population)
Lennen2015_pMA1[Lennen2015_pMA1_NaCl_row_nums, "selective_pressure"] <- "NaCl"
Lennen2015_pMA1 <- Lennen2015_pMA1 %>%
  replace(is.na(.), 0) %>%
  filter(details != "intergenic", frequency != "clone", gene != "0") %>%
  select(gene, population, selective_pressure, frequency)
Lennen2015_pMA1_out <- c(grep(",", Lennen2015_pMA1$gene), grep("/", Lennen2015_pMA1$gene))
if (length(Lennen2015_pMA1_out) > 0) {   
  Lennen2015_pMA1 <- Lennen2015_pMA1[-Lennen2015_pMA1_out,] 
  } 

write_csv(Lennen2015_pMA1, here("data_in", "for_func", "Lennen2015_pMA1.csv"))
multiple_long(paper = "Lennen2015", dataset_name = "Lennen2015_pMA1", environment = "M9", days = "5", selective_pressure = c("no NaCl", "NaCl"), 
              species = "Ecoli_K12", ploidy = "haploid", strain_info = "pMA1")


## (TL): Lennen2015_pMA7:
Lennen2015_pMA7 <- read_csv(here("data_in", "original & usable", "Lennen2015", "Lennen2015_pMA7_usable.csv"))
Lennen2015_pMA7 <- clean_names(Lennen2015_pMA7, case = "snake")
colnames(Lennen2015_pMA7) <- tolower(colnames(Lennen2015_pMA7))
### (TL): Convert to multiple_long:
Lennen2015_pMA7 <- gather(Lennen2015_pMA7, population, frequency, "m9p_ma7_1" : "na_clp_ma7_15", factor_key=TRUE)
### (TL): Add a new col, "selective_pressure". The value "0" is just a placeholder:
Lennen2015_pMA7 <- cbind(Lennen2015_pMA7, selective_pressure = 0)
### (TL): Get the rows that have M9 in the population name & save their values in "selective_pressure" as "no NaCl". Same thing for NaCl:
Lennen2015_pMA7_no_NaCl_row_nums <- grep("m9", Lennen2015_pMA7$population)
Lennen2015_pMA7[Lennen2015_pMA7_no_NaCl_row_nums, "selective_pressure"] <- "no NaCl"
Lennen2015_pMA7_NaCl_row_nums <- grep("na_cl", Lennen2015_pMA7$population)
Lennen2015_pMA7[Lennen2015_pMA7_NaCl_row_nums, "selective_pressure"] <- "NaCl"
Lennen2015_pMA7 <- Lennen2015_pMA7 %>%
  replace(is.na(.), 0) %>%
  filter(details != "intergenic", frequency != "clone", gene != "0") %>%
  select(gene, population, selective_pressure, frequency)
Lennen2015_pMA7_out <- c(grep(",", Lennen2015_pMA7$gene), grep("/", Lennen2015_pMA7$gene))
if (length(Lennen2015_pMA7_out) > 0) {   
  Lennen2015_pMA7 <- Lennen2015_pMA7[-Lennen2015_pMA7_out,] 
} 

write_csv(Lennen2015_pMA7, here("data_in", "for_func", "Lennen2015_pMA7.csv"))
multiple_long(paper = "Lennen2015", dataset_name = "Lennen2015_pMA7", environment = "M9", days = "5", selective_pressure = c("no NaCl", "NaCl"), 
              species = "Ecoli_K12", ploidy = "haploid", strain_info = "pMA7")

##############################
# Single wide:
#############################
### (TL): A template for func calling:
### (TL): single_wide(paper = "", dataset_name = "", environment = "", generations/days/flasks = "", selective_pressure = "", species = "", ploidy = "", population = c("), (strain_info = ""))

################################
# Single long:
################################
### (TL): A template for func calling:
### (TL): single_long(paper = "", dataset_name = "", environment = "", generations/days/flasks = "", selective_pressure = "", species = "", ploidy = "", (strain_info = ""))

# (Working) Long2017:
Long2017 <- read_csv(here("data_in", "original & usable", "Long2017", "Long2017_usable.csv"))
Long2017 <- clean_names(Long2017, case = "snake")
colnames(Long2017) <- tolower(colnames(Long2017))
Long2017 <- gather(Long2017, population, frequency, "ale1":"ale10", factor_key=TRUE)
Long2017 <- Long2017 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Long2017_out <- c(grep(",", Long2017$gene), grep("/", Long2017$gene))
if (length(Long2017_out) > 0) {
  Long2017 <- Long2017[-Long2017_out,] 
} 

write_csv(Long2017, here("data_in", "for_func", "Long2017.csv"))
single_long(paper = "Long2017", dataset_name = "Long2017", environment = "Glucose minimal media", days = "50", selective_pressure = "Glucose minimal media", 
            species = "Ecoli_K12", ploidy = "haploid")



# (Working) Sandberg2014:
Sandberg2014 <- read_csv(here("data_in", "original & usable", "Sandberg2014", "Sandberg2014_usable.csv"))
Sandberg2014 <- clean_names(Sandberg2014, case = "snake")
colnames(Sandberg2014) <- tolower(colnames(Sandberg2014))
Sandberg2014 <- gather(Sandberg2014, population, frequency, "ale1":"ale10", factor_key=TRUE)
Sandberg2014 <- Sandberg2014 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Sandberg2014_out <- c(grep(",", Sandberg2014$gene), grep("/", Sandberg2014$gene))
if (length(Sandberg2014_out) > 0) {
  Sandberg2014 <- Sandberg2014[-Sandberg2014_out,] 
} 

write_csv(Sandberg2014, here("data_in", "for_func", "Sandberg2014.csv"))
single_long(paper = "Sandberg2014", dataset_name = "Sandberg2014", environment = "Glucose minimal media", days = "45", 
            selective_pressure = "heat (42 degrees C)", species = "Ecoli_K12", ploidy = "haploid")



# (Working) Creamer2016:
Creamer2016 <- read_csv(here("data_in", "original & usable", "Creamer2016", "Creamer2016_usable.csv"))
Creamer2016 <- clean_names(Creamer2016, case = "snake")
colnames(Creamer2016) <- tolower(colnames(Creamer2016))
Creamer2016 <- gather(Creamer2016, population, frequency, "k0001":"k0030", factor_key=TRUE)
Creamer2016 <- Creamer2016 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Creamer2016_out <- c(grep(",", Creamer2016$gene), grep("/", Creamer2016$gene))
if (length(Creamer2016_out) > 0) {
  Creamer2016 <- Creamer2016[-Creamer2016_out,] 
} 

write_csv(Creamer2016, here("data_in", "for_func", "Creamer2016.csv"))
single_long(paper = "Creamer2016", dataset_name = "Creamer2016", environment = "LBK medium", generations = "3000", 
            selective_pressure = "benzoate", species = "Ecoli_K12", ploidy = "haploid")



# (Working) Deatherage2017:
### This paper originally had three clones sequenced for each population. Collapse the three columns and save the resulting dataset. And filter out the intergenics:
Deatherage2017 <- read_csv(here("data_in", "original & usable", "Deatherage2017", "Deatherage2017_usable.csv"))
Deatherage2017 <- clean_names(Deatherage2017, case = "snake")
colnames(Deatherage2017) <- tolower(colnames(Deatherage2017))
Deatherage2017 <- Deatherage2017 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic")
Deatherage2017 <- Deatherage2017 %>%
  transmute(gene = gene, population = population, "frequency" = `f1_i1_r1` + `f1_i2_r1` + `f1_i3_r1`)
Deatherage2017_out <- c(grep(",", Deatherage2017$gene), grep("/", Deatherage2017$gene))
if (length(Deatherage2017_out) > 0) {   
  Deatherage2017 <- Deatherage2017[-Deatherage2017_out,] 
  } 

write_csv(Deatherage2017, here("data_in", "for_func", "Deatherage2017.csv"))
single_long(paper = "Deatherage2017", dataset_name = "Deatherage2017", environment = "Glucose minimal medium", generations = "2000", 
            selective_pressure = "Glucose minimal medium", species = "Ecoli_K12", ploidy = "haploid")



# (TL) McCloskey2018:
McCloskey2018 <- read_csv(here("data_in", "original & usable", "McCloskey2018", "McCloskey2018_usable.csv"))
McCloskey2018 <- clean_names(McCloskey2018, case = "snake")
colnames(McCloskey2018) <- tolower(colnames(McCloskey2018))
McCloskey2018 <- McCloskey2018 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic", gene != "0") %>%
  select(gene, population, frequency)
McCloskey2018_out <- c(grep(",", McCloskey2018$gene), grep("/", McCloskey2018$gene))
if (length(McCloskey2018_out) > 0) {   
  McCloskey2018 <- McCloskey2018[-McCloskey2018_out,] 
}
            
## (TL) McCloskey2018_gnd:
McCloskey2018_gnd_rows <- grep("gnd", McCloskey2018$population)
McCloskey2018_gnd <- McCloskey2018[McCloskey2018_gnd_rows,]

write_csv(McCloskey2018_gnd, here("data_in", "for_func", "McCloskey2018_gnd.csv"))
single_long(paper = "McCloskey2018", dataset_name = "McCloskey2018_gnd", environment = "M9 minimal medium", days = "35", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid", strain_info = "gnd")


## (TL) McCloskey2018_pgi:
McCloskey2018_pgi_rows <- grep("pgi", McCloskey2018$population)
McCloskey2018_pgi <- McCloskey2018[McCloskey2018_pgi_rows,]

write_csv(McCloskey2018_pgi, here("data_in", "for_func", "McCloskey2018_pgi.csv"))
single_long(paper = "McCloskey2018", dataset_name = "McCloskey2018_pgi", environment = "M9 minimal medium", days = "35", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid", strain_info = "pgi")


## (TL) McCloskey2018_ptsHIcrr:
McCloskey2018_ptsHIcrr_rows <- grep("ptsHIcrr", McCloskey2018$population)
McCloskey2018_ptsHIcrr <- McCloskey2018[McCloskey2018_ptsHIcrr_rows,]

write_csv(McCloskey2018_ptsHIcrr, here("data_in", "for_func", "McCloskey2018_ptsHIcrr.csv"))
single_long(paper = "McCloskey2018", dataset_name = "McCloskey2018_ptsHIcrr", environment = "M9 minimal medium", days = "35", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid", strain_info ="ptsHIcrr")


## (TL) McCloskey2018_sdhCB:
McCloskey2018_sdhCB_rows <- grep("sdhCB", McCloskey2018$population)
McCloskey2018_sdhCB <- McCloskey2018[McCloskey2018_sdhCB_rows,]

write_csv(McCloskey2018_sdhCB, here("data_in", "for_func", "McCloskey2018_sdhCB.csv"))
single_long(paper = "McCloskey2018", dataset_name = "McCloskey2018_sdhCB", environment = "M9 minimal medium", days = "35", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid", strain_info = "sdhCB")


## (TL) McCloskey2018_tpiAE:
McCloskey2018_tpiAE_rows <- grep("tpiAE", McCloskey2018$population)
McCloskey2018_tpiAE <- McCloskey2018[McCloskey2018_tpiAE_rows,]

write_csv(McCloskey2018_tpiAE, here("data_in", "for_func", "McCloskey2018_tpiAE.csv"))
single_long(paper = "McCloskey2018", dataset_name = "McCloskey2018_tpiAE", environment = "M9 minimal medium", days = "35", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid", strain_info = "tpiAE")


## (TL) McCloskey2018_Evo:
McCloskey2018_Evo_rows <- grep("Evo", McCloskey2018$population)
McCloskey2018_Evo <- McCloskey2018[McCloskey2018_Evo_rows,]

write_csv(McCloskey2018_Evo, here("data_in", "for_func", "McCloskey2018_Evo.csv"))
single_long(paper = "McCloskey2018", dataset_name = "McCloskey2018_Evo", environment = "M9 minimal medium", days = "35", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid", strain_info = "Evo")



# (TL): Sandra2016:
Sandra2016 <- read_csv(here("data_in", "original & usable", "Sandra2016", "Sandra2016_usable.csv"))
Sandra2016 <- clean_names(Sandra2016, case = "snake")
colnames(Sandra2016) <- tolower(colnames(Sandra2016))
Sandra2016 <- Sandra2016 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Sandra2016_out <- c(grep(",", Sandra2016$gene), grep("/", Sandra2016$gene))
if (length(Sandra2016_out) > 0) {   
  Sandra2016 <- Sandra2016[-Sandra2016_out,] 
} 

write_csv(Sandra2016, here("data_in", "for_func", "Sandra2016.csv"))
single_long(paper = "Sandra2016", dataset_name = "Sandra2016", environment = "M9 minimal medium", days = "8", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid")


# (TL): Khare2015, but just with pvdJ spent media (other 2 environments only has 2 populations):
Khare2015 <- read_csv(here("data_in", "original & usable", "Khare2015", "Khare2015_usable.csv"))
Khare2015 <- clean_names(Khare2015, case = "snake")
colnames(Khare2015) <- tolower(colnames(Khare2015))
Khare2015 <- Khare2015 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
### (TL): (190703) New indicators of multiple genes involved: "[", "]", & "-".
Khare2015_out <- c(grep(",", Khare2015$gene), grep("/", Khare2015$gene))
if (length(Khare2015_out) > 0) {   
  Khare2015 <- Khare2015[-Khare2015_out,] 
} 

write_csv(Khare2015, here("data_in", "for_func", "Khare2015.csv"))
single_long(paper = "Khare2015", dataset_name = "Khare2015", environment = "LB", days = "15", 
            selective_pressure = "pvdJ spent media", species = "Ecoli_K12", ploidy = "haploid")



# (TL): Sandberg2017, subsetted by environment. 
## (TL): Sandberg2017_ac:
Sandberg2017_ac <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_ac_usable.csv"))
Sandberg2017_ac <- clean_names(Sandberg2017_ac, case = "snake")
colnames(Sandberg2017_ac) <- tolower(colnames(Sandberg2017_ac))
Sandberg2017_ac <- gather(Sandberg2017_ac, population, frequency, "a1_f56_i1_r1" : "a4_f55_i1_r1", factor_key=TRUE)
Sandberg2017_ac <- Sandberg2017_ac %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Sandberg2017_ac_out <- c(grep(",", Sandberg2017_ac$gene), grep("/", Sandberg2017_ac$gene))
if (length(Sandberg2017_ac_out) > 0) {   
  Sandberg2017_ac <- Sandberg2017_ac[-Sandberg2017_ac_out,] 
} 

write_csv(Sandberg2017_ac, here("data_in", "for_func", "Sandberg2017_ac.csv"))
single_long(paper = "Sandberg2017", dataset_name = "Sandberg2017_ac", environment = "M9 minimal medium + acetate", generations = "1000", 
            selective_pressure = "acetate", species = "Ecoli_K12", ploidy = "haploid")


## (TL): Sandberg2017_glu_ac:
Sandberg2017_glu_ac <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_ac_usable.csv"))
Sandberg2017_glu_ac <- clean_names(Sandberg2017_glu_ac, case = "snake")
colnames(Sandberg2017_glu_ac) <- tolower(colnames(Sandberg2017_glu_ac))
### (TL): transmute() here is used to collapse the clones & replicates from the same population.
Sandberg2017_glu_ac <- Sandberg2017_glu_ac %>% 
  transmute(gene = gene, details = details, 
            a7 = rowSums(Sandberg2017_glu_ac[, 4:15]), a8 = rowSums(Sandberg2017_glu_ac[, 16:28]), 
            a9 = rowSums(Sandberg2017_glu_ac[, 29:38]))
Sandberg2017_glu_ac <- gather(Sandberg2017_glu_ac, population, frequency, "a7" : "a9", factor_key=TRUE)
Sandberg2017_glu_ac <- Sandberg2017_glu_ac %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Sandberg2017_glu_ac_out <- c(grep(",", Sandberg2017_glu_ac$gene), grep("/", Sandberg2017_glu_ac$gene))
if (length(Sandberg2017_glu_ac_out) > 0) {   
  Sandberg2017_glu_ac <- Sandberg2017_glu_ac[-Sandberg2017_glu_ac_out,] 
} 

write_csv(Sandberg2017_glu_ac, here("data_in", "for_func", "Sandberg2017_glu_ac.csv"))
single_long(paper = "Sandberg2017", dataset_name = "Sandberg2017_glu_ac", environment = "M9 minimal medium + glucose + acetate", generations = "650", 
            selective_pressure = "glucose + acetate", species = "Ecoli_K12", ploidy = "haploid")


## (TL): Sandberg2017_glu_gly:
Sandberg2017_glu_gly <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_gly_usable.csv"))
Sandberg2017_glu_gly <- clean_names(Sandberg2017_glu_gly, case = "snake")
colnames(Sandberg2017_glu_gly) <- tolower(colnames(Sandberg2017_glu_gly))
Sandberg2017_glu_gly <- Sandberg2017_glu_gly %>% 
  transmute(gene = gene, details = details, 
            a4 = rowSums(Sandberg2017_glu_gly[, 4:7]), a5 = rowSums(Sandberg2017_glu_gly[, 8:12]), 
            a6 = rowSums(Sandberg2017_glu_gly[, 13:16]))
Sandberg2017_glu_gly <- gather(Sandberg2017_glu_gly, population, frequency, "a4" : "a6", factor_key=TRUE)
Sandberg2017_glu_gly <- Sandberg2017_glu_gly %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Sandberg2017_glu_gly_out <- c(grep(",", Sandberg2017_glu_gly$gene), grep("/", Sandberg2017_glu_gly$gene))
if (length(Sandberg2017_glu_gly_out) > 0) {   
  Sandberg2017_glu_gly <- Sandberg2017_glu_gly[-Sandberg2017_glu_gly_out,] 
} 

write_csv(Sandberg2017_glu_gly, here("data_in", "for_func", "Sandberg2017_glu_gly.csv"))
single_long(paper = "Sandberg2017", dataset_name = "Sandberg2017_glu_gly", environment = "M9 minimal medium + glucose + glycerol", generations = "1170", 
            selective_pressure = "glucose + glycerol", species = "Ecoli_K12", ploidy = "haploid")

## (TL): Sandberg2017_glu_xyl:
Sandberg2017_glu_xyl <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_xyl_usable.csv"))
Sandberg2017_glu_xyl <- clean_names(Sandberg2017_glu_xyl, case = "snake")
colnames(Sandberg2017_glu_xyl) <- tolower(colnames(Sandberg2017_glu_xyl))
Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl %>% 
  transmute(gene = gene, details = details, 
            a0 = rowSums(Sandberg2017_glu_xyl[, 3]), a1 = rowSums(Sandberg2017_glu_xyl[, 4:13]), 
            a2 = rowSums(Sandberg2017_glu_xyl[, 14:19]), a3 = rowSums(Sandberg2017_glu_xyl[, 20:ncol(Sandberg2017_glu_xyl)]))
Sandberg2017_glu_xyl <- gather(Sandberg2017_glu_xyl, population, frequency, "a1" : "a3", factor_key=TRUE)
Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Sandberg2017_glu_xyl_out <- c(grep(",", Sandberg2017_glu_xyl$gene), grep("/", Sandberg2017_glu_xyl$gene))
if (length(Sandberg2017_glu_xyl_out) > 0) {   
  Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl[-Sandberg2017_glu_xyl_out,] 
} 

write_csv(Sandberg2017_glu_xyl, here("data_in", "for_func", "Sandberg2017_glu_xyl.csv"))
single_long(paper = "Sandberg2017", dataset_name = "Sandberg2017_glu_xyl", environment = "M9 minimal medium + glucose + xylose", generations = "1180", 
            selective_pressure = "glucose + xylose", species = "Ecoli_K12", ploidy = "haploid")


## (TL): Sandberg2017_xyl:
Sandberg2017_xyl <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_xyl_usable.csv"))
Sandberg2017_xyl <- clean_names(Sandberg2017_xyl, case = "snake")
colnames(Sandberg2017_xyl) <- tolower(colnames(Sandberg2017_xyl))
Sandberg2017_xyl <- gather(Sandberg2017_xyl, population, frequency, "a1_f116_i1_r1" : "a4_f113_i1_r1", factor_key=TRUE)
Sandberg2017_xyl <- Sandberg2017_xyl %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Sandberg2017_xyl_out <- c(grep(",", Sandberg2017_xyl$gene), grep("/", Sandberg2017_xyl$gene))
if (length(Sandberg2017_xyl_out) > 0) {   
  Sandberg2017_xyl <- Sandberg2017_xyl[-Sandberg2017_xyl_out,] 
} 

write_csv(Sandberg2017_xyl, here("data_in", "for_func", "Sandberg2017_xyl.csv"))
single_long(paper = "Sandberg2017", dataset_name = "Sandberg2017_xyl", environment = "M9 minimal medium + xylose", generations = "1000", 
            selective_pressure = "xylose", species = "Ecoli_K12", ploidy = "haploid")



# (TL): Griffith2019, subsetted by pH:
## (TL): Griffith2019_pH_6.5:
Griffith2019_pH_6.5 <- read_csv(here("data_in", "original & usable", "Griffith2019", "Griffith2019_pH_6.5_usable.csv"))
Griffith2019_pH_6.5 <- clean_names(Griffith2019_pH_6.5, case = "snake")
colnames(Griffith2019_pH_6.5) <- tolower(colnames(Griffith2019_pH_6.5))
Griffith2019_pH_6.5 <- gather(Griffith2019_pH_6.5, population, frequency, "c_a1_1" : "c_h5_1", factor_key=TRUE)
Griffith2019_pH_6.5 <- Griffith2019_pH_6.5 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Griffith2019_pH_6.5_out <- c(grep(",", Griffith2019_pH_6.5$gene), grep("/", Griffith2019_pH_6.5$gene))
if (length(Griffith2019_pH_6.5_out) > 0) {   
  Griffith2019_pH_6.5 <- Griffith2019_pH_6.5[-Griffith2019_pH_6.5_out,] 
} 

write_csv(Griffith2019_pH_6.5, here("data_in", "for_func", "Griffith2019_pH_6.5.csv"))
single_long(paper = "Griffith2019", dataset_name = "Griffith2019_pH_6.5", environment = "Modified LBK medium", generations = "1000", 
            selective_pressure = "pH 6.5", species = "Ecoli_K12", ploidy = "haploid")


## (TL): Griffith2019_pH_8:
Griffith2019_pH_8 <- read_csv(here("data_in", "original & usable", "Griffith2019", "Griffith2019_pH_8_usable.csv"))
Griffith2019_pH_8 <- clean_names(Griffith2019_pH_8, case = "snake")
colnames(Griffith2019_pH_8) <- tolower(colnames(Griffith2019_pH_8))
Griffith2019_pH_8 <- gather(Griffith2019_pH_8, population, frequency, "c_g7_1" : "c_d11_1", factor_key=TRUE)
Griffith2019_pH_8 <- Griffith2019_pH_8 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Griffith2019_pH_8_out <- c(grep(",", Griffith2019_pH_8$gene), grep("/", Griffith2019_pH_8$gene))
if (length(Griffith2019_pH_8_out) > 0) {   
  Griffith2019_pH_8 <- Griffith2019_pH_8[-Griffith2019_pH_8_out,] 
} 

write_csv(Griffith2019_pH_8, here("data_in", "for_func", "Griffith2019_pH_8.csv"))
single_long(paper = "Griffith2019", dataset_name = "Griffith2019_pH_8", environment = "Modified LBK medium", generations = "1000", 
            selective_pressure = "pH 8", species = "Ecoli_K12", ploidy = "haploid")



# (TL): Avrani2017:
Avrani2017 <- read_csv(here("data_in", "original & usable", "Avrani2017", "Avrani2017_usable.csv"))
Avrani2017 <- clean_names(Avrani2017, case = "snake")
colnames(Avrani2017) <- tolower(colnames(Avrani2017))
Avrani2017 <- Avrani2017 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Avrani2017_out <- c(grep(",", Avrani2017$gene), grep("/", Avrani2017$gene))
if (length(Avrani2017_out) > 0) {   
  Avrani2017 <- Avrani2017[-Avrani2017_out,] 
} 

write_csv(Avrani2017, here("data_in", "for_func", "Avrani2017.csv"))
single_long(paper = "Avrani2017", dataset_name = "Avrani2017", environment = "LB", days = "11", 
            selective_pressure = "LB", species = "Ecoli_K12", ploidy = "haploid")



# (TL): Charusanti2010:
Charusanti2010 <- read_csv(here("data_in", "original & usable", "Charusanti2010", "Charusanti2010_usable.csv"))
Charusanti2010 <- clean_names(Charusanti2010, case = "snake")
colnames(Charusanti2010) <- tolower(colnames(Charusanti2010))
Charusanti2010 <- gather(Charusanti2010, population, frequency, "a1" : "a10", factor_key=TRUE)
Charusanti2010 <- Charusanti2010 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Charusanti2010_out <- c(grep(",", Charusanti2010$gene), grep("/", Charusanti2010$gene))
if (length(Charusanti2010_out) > 0) {   
  Charusanti2010 <- Charusanti2010[-Charusanti2010_out,] 
} 

write_csv(Charusanti2010, here("data_in", "for_func", "Charusanti2010.csv"))
single_long(paper = "Charusanti2010", dataset_name = "Charusanti2010", environment = "M9 minimal medium", days = "50", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid")


### (TL): 190623, 4p.
# (TL): Anand2019 (only 1 founder had 3 replicates evolving in the same condition)(gene names after underscores indicate genes missing):
Anand2019_menF_entC_ubiC <- read_csv(here("data_in", "original & usable", "Anand2019", "Anand2019_menF_entC_ubiC_usable.csv"))
Anand2019_menF_entC_ubiC <- clean_names(Anand2019_menF_entC_ubiC, case = "snake")
colnames(Anand2019_menF_entC_ubiC) <- tolower(colnames(Anand2019_menF_entC_ubiC))
Anand2019_menF_entC_ubiC <- Anand2019_menF_entC_ubiC %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Anand2019_menF_entC_ubiC_out <- c(grep(",", Anand2019_menF_entC_ubiC$gene), grep("/", Anand2019_menF_entC_ubiC$gene))
if (length(Anand2019_menF_entC_ubiC_out) > 0) {   
  Anand2019_menF_entC_ubiC <- Anand2019_menF_entC_ubiC[-Anand2019_menF_entC_ubiC_out,] 
} 

write_csv(Anand2019_menF_entC_ubiC, here("data_in", "for_func", "Anand2019_menF_entC_ubiC.csv"))
single_long(paper = "Anand2019", dataset_name = "Anand2019_menF_entC_ubiC", environment = "M9 minimal medium", days = "41", 
            selective_pressure = "iron (FeSO4)", species = "Ecoli_K12", ploidy = "haploid", strain_info = "Missing menF, entC, & ubiC genes") 



# (TL): Tenaillon2012:
Tenaillon2012 <- read_csv(here("data_in", "original & usable", "Tenaillon2012", "Tenaillon2012_usable.csv"))
Tenaillon2012 <- clean_names(Tenaillon2012, case = "snake")
colnames(Tenaillon2012) <- tolower(colnames(Tenaillon2012))
Tenaillon2012 <- gather(Tenaillon2012, population, frequency, "a0_f1_i1_r1": "a143_f1_i1_r1", factor_key=TRUE)
Tenaillon2012 <- Tenaillon2012 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Tenaillon2012_out <- c(grep(",", Tenaillon2012$gene), grep("/", Tenaillon2012$gene))
if (length(Tenaillon2012_out) > 0) {   
  Tenaillon2012 <- Tenaillon2012[-Tenaillon2012_out,] 
} 

write_csv(Tenaillon2012, here("data_in", "for_func", "Tenaillon2012.csv"))
single_long(paper = "Tenaillon2012", dataset_name = "Tenaillon2012", environment = "Davis minimal medium", generations = "2000", 
            selective_pressure = "heat (42.2 degrees C)", species = "Ecoli_K12", ploidy = "haploid")



### (TL): I TLed to run Tenaillon2012 with single_wide(), but to no avail.
### (TL): A likely reason the "undefined columns selected" error keeps popping up when trying to run with single_wide:
### (TL): The column names have spaces, which makes them syntactically incorrect.
### (TL): When running the code, something from the inner workings of the functions used to create single_wide() automatically changes column names to syntactically valid names in order to meet the requirements of some deep-level functions.
### (TL); Then, what happens is later on in single_wide(), the dataframe num_parallel TLes to find population columns whose names were entered BEFORE they were turned into syntactically correct names!
### (TL): Since the population columns now have new, syntax-abiding names, num_parallel can't find the old names, hence the "undefined columns selected" error.
### (TL): So, to tackle this problem, I use make.names(), which makes syntactically valid names from character vectors:
# names(Tenaillon2012) <- make.names(names(Tenaillon2012), unique = TRUE)
# if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Tenaillon2012, here("data_in", "for_func", "Tenaillon2012.csv"))
### (TL): The "undefined columns selected" problem is solved, but the "need finite 'ylim' values" appears at line 101 of single_wide().
# single_wide("Tenaillon2012", names(Tenaillon2012)[3:ncol(Tenaillon2012)], "Davis minimal medium", "Ecoli_K12")



# (Working) Wannier2018:
### This paper originally had two clones sequenced for each population. Collapse the two columns and save the resulting dataset. 
Wannier2018 <- read_csv(here("data_in", "original & usable", "Wannier2018", "Wannier2018_usable.csv"))
Wannier2018 <- clean_names(Wannier2018, case = "snake")
colnames(Wannier2018) <- tolower(colnames(Wannier2018))
Wannier2018 <- Wannier2018 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic")
Wannier2018 <- Wannier2018 %>%
  transmute(gene = gene, a1 = `a1_f1_i1_r1`+ `a1_f1_i2_r1`, a2 = `a2_f1_i1_r1`+ `a2_f1_i2_r1`, a3 = `a3_f1_i1_r1`+ `a3_f1_i2_r1`, 
            a4 = `a4_f1_i1_r1`+ `a4_f1_i2_r1`, a5 = `a5_f1_i1_r1`+ `a5_f1_i2_r1`, a6 = `a6_f1_i1_r1`+ `a6_f1_i2_r1`,  
            a7 = `a7_f1_i1_r1`+ `a7_f1_i2_r1`, a8 = `a8_f1_i1_r1`+ `a8_f1_i2_r1`, a9 = `a9_f1_i1_r1`+ `a9_f1_i2_r1`, 
            a10 = `a10_f1_i1_r1`+ `a10_f1_i2_r1`, a11 = `a11_f1_i1_r1`+ `a11_f1_i2_r1`, a12 = `a12_f1_i1_r1`+ `a12_f1_i2_r1`,  
            a13 = `a13_f1_i1_r1`+ `a13_f1_i2_r1`, a14 = `a14_f1_i1_r1`+ `a14_f1_i2_r1`) 
### (TL): MH's code didn't work (see below), so I converted the dataset to single long & analyze w/ single_long().
Wannier2018 <- gather(Wannier2018, population, frequency, "a1" : "a14", factor_key=TRUE)
Wannier2018 <- Wannier2018 %>%
  select(gene, population, frequency)
Wannier2018_out <- c(grep(",", Wannier2018$gene), grep("/", Wannier2018$gene))
if (length(Wannier2018_out) > 0) {   
  Wannier2018 <- Wannier2018[-Wannier2018_out,] 
} 

write_csv(Wannier2018, here("data_in", "for_func", "Wannier2018.csv"))
single_long(paper = "Wannier2018", dataset_name = "Wannier2018", environment = "M9 minimal medium", generations = "1100", 
            selective_pressure = "M9 minimal medium", species = "Ecoli_K12", ploidy = "haploid")

# (TL) MH's code for Wannier2018:
### (TL) Extract just the population columns from Wannier2018:
# Wannier2018.maTLx <- Wannier2018 %>%
#   select(population) %>%
#   as.maTLx(.)
### (TL) Replace all > 1 values with 1:
# Wannier2018.maTLx[Wannier2018.maTLx > 0] <-1
### (TL) Recombine the new population columns with genes:
# Wannier2018 <- cbind(gene = Wannier2018$gene, as.data.frame(Wannier2018.maTLx))
# if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Wannier2018, here("data_in", "Wannier2018.csv"))
### (TL) The single_wide() call isn't working - "Error in plot.window(...) : need finite 'ylim' values". 
### (TL) Error traces back to line 101 of single_wide.R, "estimate_pa(full_maTLx, ndigits = 4, show.plot = T)".
# single_wide("Wannier2018", c("A1" , "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14"), "Glucose Mininal Medium", "Ecoli_K12")



# (TL) KuzdzalFick2018:
KuzdzalFick2018 <- read_csv(here("data_in", "original & usable", "KuzdzalFick2018", "KuzdzalFick2018_usable.csv"))
KuzdzalFick2018 <- clean_names(KuzdzalFick2018, case = "snake")
colnames(KuzdzalFick2018) <- tolower(colnames(KuzdzalFick2018))
KuzdzalFick2018 <- KuzdzalFick2018 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
KuzdzalFick2018_out <- c(grep(",", KuzdzalFick2018$gene), grep("/", KuzdzalFick2018$gene))
if (length(KuzdzalFick2018_out) > 0) {   
  KuzdzalFick2018 <- KuzdzalFick2018[-KuzdzalFick2018_out,] 
} 

write_csv(KuzdzalFick2018, here("data_in", "for_func", "KuzdzalFick2018.csv"))
single_long(paper = "KuzdzalFick2018", dataset_name = "KuzdzalFick2018", environment = "YPD", days = "28", 
            selective_pressure = "M9 minimal medium", species = "Sac", ploidy = "haploid")



# (TL): Boyle2017:
Boyle2017 <- read_csv(here("data_in", "original & usable", "Boyle2017", "Boyle2017_usable.csv"))
Boyle2017 <- clean_names(Boyle2017, case = "snake")
colnames(Boyle2017) <- tolower(colnames(Boyle2017))
Boyle2017 <- Boyle2017 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Boyle2017_out <- c(grep(",", Boyle2017$gene), grep("/", Boyle2017$gene))
if (length(Boyle2017_out) > 0) {   
  Boyle2017 <- Boyle2017[-Boyle2017_out,] 
} 

write_csv(Boyle2017, here("data_in", "for_func", "Boyle2017.csv"))
single_long(paper = "Boyle2017", dataset_name = "Boyle2017", environment = "LB", generations = "49", 
            selective_pressure = "LB", species = "P_aeruginosa_PA14", ploidy = "haploid")



# (TL): Conrad2009, subset by experiment:
## (TL): Conrad2009_A_to_E:
Conrad2009_A_to_E <- read_csv(here("data_in", "original & usable", "Conrad2009", "Conrad2009_A_to_E_usable.csv"))
Conrad2009_A_to_E <- clean_names(Conrad2009_A_to_E, case = "snake")
colnames(Conrad2009_A_to_E) <- tolower(colnames(Conrad2009_A_to_E))
Conrad2009_A_to_E <- Conrad2009_A_to_E %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Conrad2009_A_to_E_out <- c(grep(",", Conrad2009_A_to_E$gene), grep("/", Conrad2009_A_to_E$gene))
if (length(Conrad2009_A_to_E_out) > 0) {   
  Conrad2009_A_to_E <- Conrad2009_A_to_E[-Conrad2009_A_to_E_out,] 
} 

write_csv(Conrad2009_A_to_E, here("data_in", "for_func", "Conrad2009_A_to_E.csv"))
single_long(paper = "Conrad2009", dataset_name = "Conrad2009_A_to_E", environment = "M9", generations = "1100", 
            selective_pressure = "Lactate", species = "Ecoli_K12", ploidy = "haploid")


## (TL): Conrad2009_F_to_K:
Conrad2009_F_to_K <- read_csv(here("data_in", "original & usable", "Conrad2009", "Conrad2009_F_to_K_usable.csv"))
Conrad2009_F_to_K <- clean_names(Conrad2009_F_to_K, case = "snake")
colnames(Conrad2009_F_to_K) <- tolower(colnames(Conrad2009_F_to_K))
Conrad2009_F_to_K <- Conrad2009_F_to_K %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Conrad2009_F_to_K_out <- c(grep(",", Conrad2009_F_to_K$gene), grep("/", Conrad2009_F_to_K$gene))
if (length(Conrad2009_F_to_K_out) > 0) {   
  Conrad2009_F_to_K <- Conrad2009_F_to_K[-Conrad2009_F_to_K_out,] 
} 

write_csv(Conrad2009_F_to_K, here("data_in", "for_func", "Conrad2009_F_to_K.csv"))
single_long(paper = "Conrad2009", dataset_name = "Conrad2009_F_to_K", environment = "M9", generations = "750", 
            selective_pressure = "Lactate", species = "Ecoli_K12", ploidy = "haploid")



# (TL): Echenique2019, subset by founder [use strain_info]:
## (TL): Echenique2019_ADE2:
Echenique2019_ADE2 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_ADE2_usable.csv"))
Echenique2019_ADE2 <- clean_names(Echenique2019_ADE2, case = "snake")
colnames(Echenique2019_ADE2) <- tolower(colnames(Echenique2019_ADE2))
Echenique2019_ADE2 <- Echenique2019_ADE2 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_ADE2_out <- c(grep(",", Echenique2019_ADE2$gene), grep("/", Echenique2019_ADE2$gene))
if (length(Echenique2019_ADE2_out) > 0) {   
  Echenique2019_ADE2 <- Echenique2019_ADE2[-Echenique2019_ADE2_out,] 
} 

write_csv(Echenique2019_ADE2, here("data_in", "for_func", "Echenique2019_ADE2.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_ADE2", strain_info = "ADE2 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_BMH1:
Echenique2019_BMH1 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_BMH1_usable.csv"))
Echenique2019_BMH1 <- clean_names(Echenique2019_BMH1, case = "snake")
colnames(Echenique2019_BMH1) <- tolower(colnames(Echenique2019_BMH1))
Echenique2019_BMH1 <- Echenique2019_BMH1 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_BMH1_out <- c(grep(",", Echenique2019_BMH1$gene), grep("/", Echenique2019_BMH1$gene))
if (length(Echenique2019_BMH1_out) > 0) {   
  Echenique2019_BMH1 <- Echenique2019_BMH1[-Echenique2019_BMH1_out,] 
} 

write_csv(Echenique2019_BMH1, here("data_in", "for_func", "Echenique2019_BMH1.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_BMH1", strain_info = "BMH1 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_COQ2:
Echenique2019_COQ2 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_COQ2_usable.csv"))
Echenique2019_COQ2 <- clean_names(Echenique2019_COQ2, case = "snake")
colnames(Echenique2019_COQ2) <- tolower(colnames(Echenique2019_COQ2))
Echenique2019_COQ2 <- Echenique2019_COQ2 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_COQ2_out <- c(grep(",", Echenique2019_COQ2$gene), grep("/", Echenique2019_COQ2$gene))
if (length(Echenique2019_COQ2_out) > 0) {   
  Echenique2019_COQ2 <- Echenique2019_COQ2[-Echenique2019_COQ2_out,] 
} 

write_csv(Echenique2019_COQ2, here("data_in", "for_func", "Echenique2019_COQ2.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_COQ2", strain_info = "COQ2 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_COX6:
Echenique2019_COX6 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_COX6_usable.csv"))
Echenique2019_COX6 <- clean_names(Echenique2019_COX6, case = "snake")
colnames(Echenique2019_COX6) <- tolower(colnames(Echenique2019_COX6))
Echenique2019_COX6 <- Echenique2019_COX6 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_COX6_out <- c(grep(",", Echenique2019_COX6$gene), grep("/", Echenique2019_COX6$gene))
if (length(Echenique2019_COX6_out) > 0) {   
  Echenique2019_COX6 <- Echenique2019_COX6[-Echenique2019_COX6_out,] 
} 

write_csv(Echenique2019_COX6, here("data_in", "for_func", "Echenique2019_COX6.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_COX6", strain_info = "COX6 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_CTF19:
Echenique2019_CTF19 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_CTF19_usable.csv"))
Echenique2019_CTF19 <- clean_names(Echenique2019_CTF19, case = "snake")
colnames(Echenique2019_CTF19) <- tolower(colnames(Echenique2019_CTF19))
Echenique2019_CTF19 <- Echenique2019_CTF19 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_CTF19_out <- c(grep(",", Echenique2019_CTF19$gene), grep("/", Echenique2019_CTF19$gene))
if (length(Echenique2019_CTF19_out) > 0) {   
  Echenique2019_CTF19 <- Echenique2019_CTF19[-Echenique2019_CTF19_out,] 
} 

write_csv(Echenique2019_CTF19, here("data_in", "for_func", "Echenique2019_CTF19.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_CTF19", strain_info = "CTF19 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_ELP4:
Echenique2019_ELP4 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_ELP4_usable.csv"))
Echenique2019_ELP4 <- clean_names(Echenique2019_ELP4, case = "snake")
colnames(Echenique2019_ELP4) <- tolower(colnames(Echenique2019_ELP4))
Echenique2019_ELP4 <- Echenique2019_ELP4 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_ELP4_out <- c(grep(",", Echenique2019_ELP4$gene), grep("/", Echenique2019_ELP4$gene))
if (length(Echenique2019_ELP4_out) > 0) {   
  Echenique2019_ELP4 <- Echenique2019_ELP4[-Echenique2019_ELP4_out,] 
} 

write_csv(Echenique2019_ELP4, here("data_in", "for_func", "Echenique2019_ELP4.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_ELP4", strain_info = "ELP4 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_HXK2:
Echenique2019_HXK2 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_HXK2_usable.csv"))
Echenique2019_HXK2 <- clean_names(Echenique2019_HXK2, case = "snake")
colnames(Echenique2019_HXK2) <- tolower(colnames(Echenique2019_HXK2))
Echenique2019_HXK2 <- Echenique2019_HXK2 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_HXK2_out <- c(grep(",", Echenique2019_HXK2$gene), grep("/", Echenique2019_HXK2$gene))
if (length(Echenique2019_HXK2_out) > 0) {   
  Echenique2019_HXK2 <- Echenique2019_HXK2[-Echenique2019_HXK2_out,] 
} 

write_csv(Echenique2019_HXK2, here("data_in", "for_func", "Echenique2019_HXK2.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_HXK2", strain_info = "HXK2 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_IML3:
Echenique2019_IML3 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_IML3_usable.csv"))
Echenique2019_IML3 <- clean_names(Echenique2019_IML3, case = "snake")
colnames(Echenique2019_IML3) <- tolower(colnames(Echenique2019_IML3))
Echenique2019_IML3 <- Echenique2019_IML3 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_IML3_out <- c(grep(",", Echenique2019_IML3$gene), grep("/", Echenique2019_IML3$gene))
if (length(Echenique2019_IML3_out) > 0) {   
  Echenique2019_IML3 <- Echenique2019_IML3[-Echenique2019_IML3_out,] 
} 

write_csv(Echenique2019_IML3, here("data_in", "for_func", "Echenique2019_IML3.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_IML3", strain_info = "IML3 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_KTI12:
Echenique2019_KTI12 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_KTI12_usable.csv"))
Echenique2019_KTI12 <- clean_names(Echenique2019_KTI12, case = "snake")
colnames(Echenique2019_KTI12) <- tolower(colnames(Echenique2019_KTI12))
Echenique2019_KTI12 <- Echenique2019_KTI12 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_KTI12_out <- c(grep(",", Echenique2019_KTI12$gene), grep("/", Echenique2019_KTI12$gene))
if (length(Echenique2019_KTI12_out) > 0) {   
  Echenique2019_KTI12 <- Echenique2019_KTI12[-Echenique2019_KTI12_out,] 
} 

write_csv(Echenique2019_KTI12, here("data_in", "for_func", "Echenique2019_KTI12.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_KTI12", strain_info = "KTI12 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_SOK2:
Echenique2019_SOK2 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_SOK2_usable.csv"))
Echenique2019_SOK2 <- clean_names(Echenique2019_SOK2, case = "snake")
colnames(Echenique2019_SOK2) <- tolower(colnames(Echenique2019_SOK2))
Echenique2019_SOK2 <- Echenique2019_SOK2 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_SOK2_out <- c(grep(",", Echenique2019_SOK2$gene), grep("/", Echenique2019_SOK2$gene))
if (length(Echenique2019_SOK2_out) > 0) {   
  Echenique2019_SOK2 <- Echenique2019_SOK2[-Echenique2019_SOK2_out,] 
} 

write_csv(Echenique2019_SOK2, here("data_in", "for_func", "Echenique2019_SOK2.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_SOK2", strain_info = "SOK2 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")


# (TL): Echenique2019_VPS29:
Echenique2019_VPS29 <- read_csv(here("data_in", "original & usable", "Echenique2019", "Echenique2019_VPS29_usable.csv"))
Echenique2019_VPS29 <- clean_names(Echenique2019_VPS29, case = "snake")
colnames(Echenique2019_VPS29) <- tolower(colnames(Echenique2019_VPS29))
Echenique2019_VPS29 <- Echenique2019_VPS29 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>%
  select(gene, population, frequency)
Echenique2019_VPS29_out <- c(grep(",", Echenique2019_VPS29$gene), grep("/", Echenique2019_VPS29$gene))
if (length(Echenique2019_VPS29_out) > 0) {   
  Echenique2019_VPS29 <- Echenique2019_VPS29[-Echenique2019_VPS29_out,] 
} 

write_csv(Echenique2019_VPS29, here("data_in", "for_func", "Echenique2019_VPS29.csv"))
single_long(paper = "Echenique2019", dataset_name = "Echenique2019_VPS29", strain_info = "VPS29 deleted", environment = "YPD", generations = "500", 
            selective_pressure = "YPD", species = "Sac", ploidy = "haploid")

####################################

# (TL) Combine all the analysis files:
file_list <- list.files("data_out/analyses", full.names = TRUE)
### (TL) Load the library "plyr" for this particular line of code only - calling this library earlier would lead to "unused arguments" error when using here().
### (TL): I think it's because the function here() is also included in the library "plyr", but it works differently, so the 2 here()'s clash.
master_analyses <- plyr::ldply(file_list, read_csv)
write_csv(master_analyses, "data_out/master_analyses.csv")

####################################
# (TL) In development: Integrate gather() to single_wide() - easier to process, less repetitive codes.
# (TL) In development: Argument "collapseMutations()" - treat all mutations in a gene as 1.
# (TL) In development: (?) Divide data_out into sub-folders, based on func used. 
