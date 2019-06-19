install.packages("Hmisc")
install.packages("plyr")
library(here)
library(dgconstraint)
library(tidyverse)
library(R.utils)
library(janitor)
library(tidyverse)
library(readr)
library(devtools)
library(Hmisc)
source("R/dgconstraint/functions/multiple_wide.R")
source("R/dgconstraint/functions/multiple_long.R")
source("R/dgconstraint/functions/single_wide.R")
source("R/dgconstraint/functions/single_long.R")

############################
# (Tri): Note: All files imported for the analysis must have "_usable.csv" at the end, for uniformity.
############################

############################
# Multiple wide:
############################

# (Working) Tenaillon2016:
### This paper originally had two clones sequenced for each timepoint. Collapse the two columns and save the resulting dataset.
Tenaillon2016 <- read_csv(here("data_in", "original & usable", "Tenaillon2016", "Tenaillon2016_usable.csv"))
Tenaillon2016 <- clean_names(Tenaillon2016)
Tenaillon2016 <- Tenaillon2016 %>% 
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") %>% 
  transmute(gene = gene, population = population, "500" = `x500_i1_r1`+`x500_i2_r1`, "1000" =`x1000_i1_r1`+`x1000_i2_r1`, "1500" =`x1500_i1_r1`+`x1500_i2_r1`,  "2000"= `x2000_i1_r1`+`x2000_i2_r1`, "5000"=`x5000_i1_r1`+`x5000_i2_r1`, "10000"= `x10000_i1_r1`+`x10000_i2_r1`, "15000"=`x15000_i1_r1`+`x15000_i2_r1`,"20000"=`x20000_i1_r1`+`x20000_i2_r1`,"30000"= `x30000_i1_r1`+`x30000_i2_r1`,"40000"=`x40000_i1_r1`+`x40000_i2_r1`,"50000"=`x50000_i1_r1`+`x50000_i2_r1`)
### (Tri): Create a new variable to store the index of rows that contain intergenic mutations not notated by breseq.
### (Tri): If there are any rows of this kind, remove them from the dataset. The if() was used because the dataset will be empty if there are no badly-notated rows (dataset(-0,) is empty!)
### (Tri): grep() is used for all hitherto known patterns indicating intergenicness in the colum "gene". Future updates could include more patterns, depending on how creative people could get with their notations.
Tenaillon2016_out <- c(grep(",", Tenaillon2016$gene), grep("/", Tenaillon2016$gene))
if (length(Tenaillon2016_out) > 0) {
  Tenaillon2016 <- Tenaillon2016[-Tenaillon2016_out,]
}

write_csv(Tenaillon2016, here("data_in", "for_func", "Tenaillon2016.csv"))
multiple_wide("Tenaillon2016", "Davis minimal medium", c("500", "1000", "1500", "2000", '5000', '10000', '20000', '30000', '40000','50000'), "Davis minimal medium", "Ecoli_K12", "haploid")



# (Working) Lang2013:
Lang2013 <- read_csv(here("data_in", "original & usable", "Lang2013", "Lang2013_usable.csv"))
Lang2013 <- Lang2013 %>%
  ### Filter out intergenics before running the code:
  filter(Gene != "Intergenic") %>% 
  rename(gene = Gene, population = Population)
Lang2013_out <- c(grep(",", Lang2013$gene), grep("/", Lang2013$gene))
if (length(Lang2013_out) > 0) {
Lang2013 <- Lang2013[-Lang2013_out,]
}

write_csv(Lang2013, here("data_in", "for_func", "Lang2013.csv"))
multiple_wide("Lang2013", "YPD", c("0","140","240","335",'415','505','585','665','745','825','910','1000'), "YPD", "Sac", "haploid")



# (Working) Sherlock2013:
Sherlock2013 <- read_csv(here("data_in", "original & usable", "Sherlock2013", "Sherlock2013_usable.csv"))
Sherlock2013 <- Sherlock2013 %>%
  filter(Gene != "Intergenic") %>% 
  rename(gene = Gene, population = Population)
# (Tri): Get rid of the "G" in generation numbers to make the "generation" column numeral:
colnames(Sherlock2013) <- gsub("G", "", colnames(Sherlock2013))
Sherlock2013_out <- c(grep(",", Sherlock2013$gene), grep("/", Sherlock2013$gene))
if (length(Sherlock2013_out) > 0) {
  Sherlock2013 <- Sherlock2013[-Sherlock2013_out,]
}

write_csv(Sherlock2013, here("data_in", "for_func", "Sherlock2013.csv"))
multiple_wide("Sherlock2013", "YPD", c("7","70", "133","196","266", "322","385","448"), "YPD", "Sac", "haploid")



# (Working) Sherlock2019:
Sherlock2019 <- read_csv(here("data_in", "original & usable", "Sherlock2019", "Sherlock2019_usable.csv"))
Sherlock2019 <- Sherlock2019 %>%
  filter(Gene != "Intergenic") %>% 
  rename(gene = Gene, population = Population)
Sherlock2019_out <- c(grep(",", Sherlock2019$gene), grep("/", Sherlock2019$gene))
if (length(Sherlock2019_out) > 0) {   
  Sherlock2019 <- Sherlock2019[-Sherlock2019_out,] 
} 
write_csv(Sherlock2019, here("data_in", "for_func", "Sherlock2019.csv"))
multiple_wide("Sherlock2019", "Davis minimal medium", c("0","50","100","150","200" ,"250","300","350","400","450","500"), "Davis minimal medium", "Ecoli_K12")



# (Working) Tonoyan2019:
Tonoyan2019 <- read_csv(here("data_in", "original & usable", "Tonoyan2019", "Tonoyan2019_usable.csv"))
Tonoyan2019 <- Tonoyan2019 %>%
  filter(Annotation != "intergenic") %>% 
  rename(gene = Gene, population = Population)

Tonoyan2019_out <- c(grep(",", Tonoyan2019$gene), grep("/", Tonoyan2019$gene))
if (length(Tonoyan2019_out) > 0) {
  Tonoyan2019 <- Tonoyan2019[-Tonoyan2019_out,] 
} 

write_csv(Tonoyan2019, here("data_in", "for_func", "Tonoyan2019.csv"))
multiple_wide("Tonoyan2019", "LB medium", c("1", "14", "20"), "LB medium", "Ecoli_ATCC")



# (Working) Wielgoss2016, subsetted by presence of mucoid:
## (Working) Wielgoss2016_mucoid:
Wielgoss2016_mucoid <- read_csv(here("data_in", "original & usable", "Wielgoss2016", "Wielgoss2016_mucoid_usable.csv"))
Wielgoss2016_mucoid <- Wielgoss2016_mucoid %>%
  filter(Gene != "Intergenic") %>% 
  rename(gene = Gene, population = Population)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Wielgoss2016_mucoid, here("data_in", "for_func", "Wielgoss2016_mucoid.csv"))
multiple_wide("Wielgoss2016_mucoid", "Lysogeny broth", c("0", "10"), "Lysogeny broth", "Ecoli_K12")


## (Working) Wielgoss2016_nonMucoid:
Wielgoss2016_nonMucoid <- read_csv(here("data_in", "original & usable", "Wielgoss2016", "Wielgoss2016_nonMucoid_usable.csv"))
Wielgoss2016_nonMucoid <- Wielgoss2016_nonMucoid %>%
  filter(Gene != "Intergenic") %>% 
  rename(gene = Gene, population = Population)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Wielgoss2016_nonMucoid, here("data_in", "for_func", "Wielgoss2016_nonMucoid.csv"))
multiple_wide("Wielgoss2016_nonMucoid", "Lysogeny broth", c("0","2", "10"), "Lysogeny broth", "Ecoli_K12")



# (Working) Sandberg2016:
Sandberg2016 <- read_csv(here("data_in", "original & usable", "Sandberg2016", "Sandberg2016_usable.csv"))
Sandberg2016 <- Sandberg2016 %>%
  filter(Gene != "Intergenic") %>% 
  rename(gene = Gene, population = Population)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Sandberg2016, here("data_in", "for_func", "Sandberg2016.csv"))
multiple_wide("Sandberg2016", "Davis minimal medium", c("Flask 23", "Flask 58", "Flask 133"), "Davis minimal medium", "Ecoli_K12")



# (Tri): Kacar2017:
Kacar2017 <- read_csv(here("data_in", "original & usable", "Kacar2017", "Kacar2017_usable.csv"))
Kacar2017 <- Kacar2017 %>% 
  filter(`Details` != "Intergenic") %>%
  rename(gene = Gene, population = Population)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Kacar2017, here("data_in", "for_func", "Kacar2017.csv"))
multiple_wide("Kacar2017", "Minimal glucose medium", c("500", "1000", "1500", "2000"), "Minimal glucose medium", "Ecoli_K12")



# (Tri): Morgenthaler2019:
Morgenthaler2019 <- read_csv(here("data_in", "original & usable", "Morgenthaler2019", "Morgenthaler2019_usable.csv"))
Morgenthaler2019 <- Morgenthaler2019 %>% 
  filter(mutation != "intergenic")
Morgenthaler2019_out <- grep("/", Morgenthaler2019$gene)
Morgenthaler2019 <- Morgenthaler2019[-Morgenthaler2019_out,]

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Morgenthaler2019, here("data_in", "for_func", "Morgenthaler2019.csv"))
multiple_wide("Morgenthaler2019", "M9", c("Day 42", "Day 50"), "M9", "Ecoli_K12")

#####################################
# Multiple long:
#####################################

# (Working) Hong2014:
Hong2014 <- read_csv(here("data_in", "original & usable", "Hong2014", "Hong2014_usable.csv"))
Hong2014 <- Hong2014 %>% 
  ### Filter out unwanted mutations:
  filter(type != 'DIP', type !='WCA', type != 'IND', type != 'AMP') %>% 
  rename(gene = Gene, population = Population)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Hong2014, here("data_in", "for_func", "Hong2014.csv"))
multiple_long("Hong2014", c("Ammonium", "Arginine", "Urea", "Allantoin"), "250", c("Ammonium", "Arginine", "Urea", "Allantoin"), "Sac")



# (Working) Jerison2017:
### For this data set there were multiple different ways to run an analysis, so each of the following code chunks will result in a different analysis.
Jerison2017 <- read_csv(here("data_in", "original & usable", "Jerison2017", "Jerison2017_usable.csv"))

### Jerison2017, subsetted by selective_pressure:
## Jerison2017_OT:
Jerison2017_OT <- Jerison2017 %>% 
  filter(Gene != "Intergenic") %>% 
  filter(selective_pressure == "OT") %>% 
  rename(gene = Gene, population = Population)
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Jerison2017_OT, here("data_in", "for_func", "Jerison2017_OT.csv"))


## Jerison2017_HT:
Jerison2017_HT <- Jerison2017 %>% 
  filter(Gene != "Intergenic") %>% 
  filter(selective_pressure == "HT") %>% 
  rename(gene = Gene, population = Population)
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Jerison2017_HT, here("data_in", "for_func", "Jerison2017_HT.csv"))


## Analyze the entire thing:
Jerison2017<-Jerison2017 %>% 
  filter(Gene != "Intergenic") %>% 
  rename(gene = Gene, population = Population)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Jerison2017, here("data_in", "for_func", "Jerison2017.csv"))
multiple_long("Jerison2017", c("YPD", "SC"), "500", c("OT", "HT"), "Sac")
single_long("Jerison2017_OT","YPD", "500", "OT", "Sac")
single_long("Jerison2017_HT","SC", "500", "HT", "Sac")



# (Working) Payen2016:
### There are a few ways to run this analysis. Whichever way you choose, you must filter out the values whose frequncies are "clone" or "frequency".
## Analyze the entire dataset, regardless of ploidy:
Payen2016 <- read_csv(here("data_in", "original & usable", "Payen2016", "Payen2016_usable.csv"))
Payen2016 <- Payen2016 %>% 
  rename(gene = Gene, population = Population) %>%
  filter(frequency != "clone") %>% 
  drop_na(frequency) 
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Payen2016, here("data_in", "for_func", "Payen2016.csv"))


## Analyze between the selective pressures for only the diploid genes:
Payen2016_dip <- Payen2016 %>% 
  filter(Ploidy == "diploid", frequency != "clone") %>% 
  drop_na(frequency) 
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Payen2016_dip, here("data_in", "for_func", "Payen2016_dip.csv"))


## Analyze between the selective pressures for only the haploid genes:
Payen2016_hap <- Payen2016 %>% 
  filter(Ploidy == "haploid", frequency != "clone") %>% 
  drop_na(frequency) 
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Payen2016_hap, here("data_in", "for_func", "Payen2016_hap.csv"))

multiple_long("Payen2016", c("phosphate", "sulfate", "glucose"), "20",  c("phosphate", "sulfate", "glucose"), "Sac")
multiple_long("Payen2016_hap", c("phosphate"), "20", c("phosphate"), "Sac")
multiple_long("Payen2016_dip", c("phosphate", "sulfate", "glucose"), "20", c("phosphate", "sulfate", "glucose"), "Sac")



# (Tri): (Working) Lennen2015, subsetted by founder (founders have different plasmids):
## (Tri): Lennen2015_pMA1:
Lennen2015_pMA1 <- read_csv(here("data_in", "original & usable", "Lennen2015", "Lennen2015_pMA1_usable.csv"))
Lennen2015_pMA1 <- Lennen2015_pMA1 %>%
  replace(is.na(.), 0)
### (Tri): Convert to multiple_long:
Lennen2015_pMA1 <- gather(Lennen2015_pMA1, population, frequency, "M9pMA1-1" : "NaClpMA1-19", factor_key=TRUE)
### (Tri): Add a column for environment. The "0" here is completely random - it's just a placeholder:
Lennen2015_pMA1 <- cbind(Lennen2015_pMA1, environment = 0)
### (Tri): Get the rows that have M9 in the population name & make that the value in "selective_pressure" (for function purposes). Same thing for NaCl:
Lennen2015_pMA1_M9_row_nums <- grep("M9", Lennen2015_pMA1$population)
Lennen2015_pMA1$selective_pressure[Lennen2015_pMA1_M9_row_nums] <- "M9"
Lennen2015_pMA1_NaCl_row_nums <- grep("NaCl", Lennen2015_pMA1$population)
Lennen2015_pMA1$selective_pressure[Lennen2015_pMA1_NaCl_row_nums] <- "NaCl"

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Lennen2015_pMA1, here("data_in", "for_func", "Lennen2015_pMA1.csv"))
multiple_long("Lennen2015_pMA1", c("M9", "NaCl"), "5 days", c("M9", "NaCl"), "Ecoli_K12")


## (Tri): Lennen2015_pMA7:
Lennen2015_pMA7 <- read_csv(here("data_in", "original & usable", "Lennen2015", "Lennen2015_pMA7_usable.csv"))
Lennen2015_pMA7 <- Lennen2015_pMA7 %>%
  replace(is.na(.), 0)
Lennen2015_pMA7 <- gather(Lennen2015_pMA7, population, frequency, "M9pMA7-1" : "NaClpMA7-15", factor_key=TRUE)
Lennen2015_pMA7 <- cbind(Lennen2015_pMA7, environment = 0)
Lennen2015_pMA7_M9_row_nums <- grep("M9", Lennen2015_pMA7$population)
Lennen2015_pMA7$selective_pressure[Lennen2015_pMA7_M9_row_nums] <- "M9"
Lennen2015_pMA7_NaCl_row_nums <- grep("NaCl", Lennen2015_pMA7$population)
Lennen2015_pMA7$selective_pressure[Lennen2015_pMA7_NaCl_row_nums] <- "NaCl"

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Lennen2015_pMA7, here("data_in", "for_func", "Lennen2015_pMA7.csv"))
multiple_long("Lennen2015_pMA7", c("M9", "NaCl"), "5 days", c("M9", "NaCl"), "Ecoli_K12")

##############################
# Single wide:
#############################

# (Working) Long2017:
Long2017 <- read_csv(here("data_in", "original & usable", "Long2017", "Long2017_usable.csv"))

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Long2017, here("data_in", "for_func", "Long2017.csv"))
single_wide("Long2017", "Glucose minimal media", "50 days", "Glucose minimal media", "Ecoli_K12", c("ALE1", "ALE2", "ALE3", "ALE4", "ALE5", "ALE6", "ALE7", "ALE8", "ALE9", "ALE10"))



# (Working) Sandberg2014:
Sandberg2014 <- read_csv(here("data_in", "original & usable", "Sandberg2014", "Sandberg2014_usable.csv"))
Sandberg2014 <- Sandberg2014 %>% 
  rename(gene = Gene)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Sandberg2014, here("data_in", "for_func", "Sandberg2014.csv"))
single_wide("Sandberg2014", "Glucose minimal media", "45 days", "Glucose minimal media", "Ecoli_K12", c("ALE1", "ALE2", "ALE3", "ALE4", "ALE5", "ALE6", "ALE7", "ALE8", "ALE9", "ALE10"))



# (Working) Creamer2016:
Creamer2016 <- read_csv(here("data_in", "original & usable", "Creamer2016", "Creamer2016_usable.csv"))
Creamer2016 <-Creamer2016 %>% 
  filter(Annotation!= "intergenic") %>% 
  remove_empty("cols") %>% 
  remove_empty("rows") %>% 
  rename(gene = Gene)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Creamer2016, here("data_in", "for_func", "Creamer2016.csv"))
single_wide("Creamer2016", "LBK medium", "2000", "LBK medium", "Ecoli_K12", c("K0001","K0011","K0002","K0003","K0014","K0015","KB026","K0022","K0023","K0006","K0007","K0010", "K0020","K0019","K0031","K0030"))

################################
# Single long:
################################

# (Working) Deatherage2017:
### This paper originally had three clones sequenced for each population. Collapse the three columns and save the resulting dataset. And filter out the intergenics:
Deatherage2017 <- read_csv(here("data_in", "original & usable", "Deatherage2017", "Deatherage2017_usable.csv"))
Deatherage2017<- Deatherage2017 %>% filter(Details != "intergenic") %>%
  replace(is.na(.), 0) %>% 
  transmute(gene = Gene, population = Population, frequency = `F1 I1 R1` + `F1 I2 R1`+`F1 I3 R1`)
   
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Deatherage2017, here("data_in", "for_func", "Deatherage2017.csv"))
single_long("Deatherage2017", "Glucose minimal medium", "2000", "Glucose minimal medium", "Ecoli_K12")



# (Working) McCloskey2018:
### Trying to find a better way with a single function to do all of the following functions:
McCloskey2018 <- read_csv(here("data_in", "original & usable", "McCloskey2018", "McCloskey2018_usable.csv"))

## For gnd, filter out everything not gnd:
McCloskey2018_gnd <- McCloskey2018 %>% 
  filter(population!="Evo04pgiEvo02EP",population!="Evo04Evo01EP",population!="Evo04Evo02EP",population!="Evo04pgiEvo01EP",population!="Evo04pgiEvo03EP",population!="Evo04pgiEvo04EP",population!="Evo04pgiEvo05EP",population!="Evo04pgiEvo06EP",population!="Evo04pgiEvo07EP",population!="Evo04pgiEvo08EP",population!="Evo04ptsHIcrrEvo01EP",population!="Evo04ptsHIcrrEvo02EP",population!="Evo04ptsHIcrrEvo03EP", population!="Evo04ptsHIcrrEvo04EP",population!= "Evo04sdhCBEvo02EP",population!="Evo04sdhCBEvo03EP",population!="Evo04tpiAEvo01EP",population!="Evo04tpiAEvo02EP",population!="Evo04tpiAEvo04EP",population!="Evo04tpiAEvo03EP",population!="Evo04sdhCBEvo01EP") %>% 
  filter(mutation_locations != "intergenic")

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(McCloskey2018_gnd, here("data_in", "for_func", "McCloskey2018_gnd.csv"))


## For pgi, filter out everything not pgi:
McCloskey2018_pgi <- McCloskey2018 %>% 
  filter(population!="Evo04Evo01EP",population!="Evo04Evo02EP",population!="Evo04gndEvo01EP",population!="Evo04gndEvo02EP",population!="Evo04gndEvo03EP",population!="Evo04ptsHIcrrEvo01EP",population!="Evo04ptsHIcrrEvo02EP",population!="Evo04ptsHIcrrEvo03EP", population!="Evo04ptsHIcrrEvo04EP",population!= "Evo04sdhCBEvo02EP",population!="Evo04sdhCBEvo03EP",population!="Evo04tpiAEvo01EP",population!="Evo04tpiAEvo02EP",population!="Evo04tpiAEvo04EP",population!="Evo04tpiAEvo03EP",population!="Evo04sdhCBEvo01EP")%>%
  filter(mutation_locations != "intergenic")

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(McCloskey2018_pgi, here("data_in", "for_func", "McCloskey2018_pgi.csv"))


## For ptsHIcrr, filter out everything not ptsHIcrr:
McCloskey2018_ptsHIcrr <- McCloskey2018 %>% 
  filter(population!="Evo04pgiEvo02EP",population!="Evo04Evo01EP",population!="Evo04Evo02EP",population!="Evo04gndEvo01EP",population!="Evo04gndEvo02EP",population!="Evo04gndEvo03EP",population!="Evo04pgiEvo01EP",population!="Evo04pgiEvo03EP",population!="Evo04pgiEvo04EP",population!="Evo04pgiEvo05EP",population!="Evo04pgiEvo06EP",population!="Evo04pgiEvo07EP",population!="Evo04pgiEvo08EP",population!= "Evo04sdhCBEvo02EP",population!="Evo04sdhCBEvo03EP",population!="Evo04tpiAEvo01EP",population!="Evo04tpiAEvo02EP",population!="Evo04tpiAEvo04EP",population!="Evo04tpiAEvo03EP",population!="Evo04sdhCBEvo01EP")%>% 
  filter(mutation_locations != "intergenic")

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(McCloskey2018_ptsHIcrr, here("data_in", "for_func", "McCloskey2018_ptsHIcrr.csv"))


## For sdhCB, filter out everything not sdhCB:
McCloskey2018_sdhCB <- McCloskey2018 %>% 
  filter(population!="Evo04pgiEvo02EP",population!="Evo04Evo01EP",population!="Evo04Evo02EP",population!="Evo04gndEvo01EP",population!="Evo04gndEvo02EP",population!="Evo04gndEvo03EP",population!="Evo04pgiEvo01EP",population!="Evo04pgiEvo03EP",population!="Evo04pgiEvo04EP",population!="Evo04pgiEvo05EP",population!="Evo04pgiEvo06EP",population!="Evo04pgiEvo07EP",population!="Evo04pgiEvo08EP",population!="Evo04ptsHIcrrEvo01EP",population!="Evo04ptsHIcrrEvo02EP",population!="Evo04ptsHIcrrEvo03EP", population!="Evo04ptsHIcrrEvo04EP",population!="Evo04tpiAEvo01EP",population!="Evo04tpiAEvo02EP",population!="Evo04tpiAEvo04EP",population!="Evo04tpiAEvo03EP")%>% 
  filter(mutation_locations != "intergenic")

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(McCloskey2018_sdhCB, here("data_in", "for_func", "McCloskey2018_sdhCB.csv"))


## For tpiAE, filter out everything not tpiAE:
McCloskey2018_tpiAE <- McCloskey2018 %>% 
  filter(population!="Evo04pgiEvo02EP",population!="Evo04Evo01EP",population!="Evo04Evo02EP",population!="Evo04gndEvo01EP",population!="Evo04gndEvo02EP",population!="Evo04gndEvo03EP",population!="Evo04pgiEvo01EP",population!="Evo04pgiEvo03EP",population!="Evo04pgiEvo04EP",population!="Evo04pgiEvo05EP",population!="Evo04pgiEvo06EP",population!="Evo04pgiEvo07EP",population!="Evo04pgiEvo08EP",population!="Evo04ptsHIcrrEvo01EP",population!="Evo04ptsHIcrrEvo02EP",population!="Evo04ptsHIcrrEvo03EP", population!="Evo04ptsHIcrrEvo04EP",population!= "Evo04sdhCBEvo02EP",population!="Evo04sdhCBEvo03EP",population!="Evo04sdhCBEvo01EP")%>% 
  filter(mutation_locations != "intergenic")

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(McCloskey2018_tpiAE, here("data_in", "for_func", "McCloskey2018_tpiAE.csv"))


## For the reference (evo), filter out everything not evo:
McCloskey2018_evo <- McCloskey2018 %>% 
  filter(population!="Evo04pgiEvo02EP",population!="Evo04gndEvo01EP",population!="Evo04gndEvo02EP",population!="Evo04gndEvo03EP",population!="Evo04pgiEvo01EP",population!="Evo04pgiEvo03EP",population!="Evo04pgiEvo04EP",population!="Evo04pgiEvo05EP",population!="Evo04pgiEvo06EP",population!="Evo04pgiEvo07EP",population!="Evo04pgiEvo08EP",population!="Evo04ptsHIcrrEvo01EP",population!="Evo04ptsHIcrrEvo02EP",population!="Evo04ptsHIcrrEvo03EP", population!="Evo04ptsHIcrrEvo04EP",population!= "Evo04sdhCBEvo02EP",population!="Evo04sdhCBEvo03EP",population!="Evo04tpiAEvo01EP",population!="Evo04tpiAEvo02EP",population!="Evo04tpiAEvo04EP",population!="Evo04tpiAEvo03EP",population!="Evo04sdhCBEvo01EP")%>%
           filter(mutation_locations != "intergenic")

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(McCloskey2018_evo, here("data_in", "for_func", "McCloskey2018_evo.csv"))


## Or just run as one function, each being an individual population:
McCloskey2018 <- McCloskey2018 %>% 
  filter(mutation_locations != "intergenic")
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(McCloskey2018, here("data_in", "for_func", "McCloskey2018.csv"))

single_long("McCloskey2018_gnd", "M9 minimal medium", "35 days" ,"M9 minimal medium", "Ecoli_K12")
single_long("McCloskey2018_pgi", "M9 minimal medium", "35 days" ,"M9 minimal medium", "Ecoli_K12")
single_long("McCloskey2018_ptsHIcrr", "M9 minimal medium", "35 days" ,"M9 minimal medium", "Ecoli_K12")
single_long("McCloskey2018_sdhCB", "M9 minimal medium", "35 days" ,"M9 minimal medium", "Ecoli_K12")
single_long("McCloskey2018_tpiAE", "M9 minimal medium", "35 days" ,"M9 minimal medium", "Ecoli_K12")
single_long("McCloskey2018_evo", "M9 minimal medium", "35 days" ,"M9 minimal medium", "Ecoli_K12")
single_long("McCloskey2018", "M9 minimal medium", "35 days" ,"M9 minimal medium", "Ecoli_K12")



# (Tri): Sandra2016:
Sandra2016 <- read_csv(here("data_in", "original & usable", "Sandra2016", "Sandra2016_usable.csv"))
Sandra2016 <- Sandra2016 %>% 
  filter(annotation != "intergenic")
### (Tri): Since this dataset doesn't have frequencies, we will assume all frequencies to be 1 for the time being (manually added to Excel):
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Sandra2016, here("data_in", "for_func", "Sandra2016.csv"))
single_long("Sandra2016", "M9 minimal medium", "8 days", "M9 minimal medium", "Ecoli_K12")



# (Tri): Khare2015, but just with pvdJ spent media (other 2 environments only has 2 populations):
Khare2015 <- read_csv(here("data_in", "original & usable", "Khare2015", "Khare2015_usable.csv"))

Khare2015 <- Khare2015 %>% 
  filter(annotation != "intergenic")
## (Tri): Since this dataset doesn't have frequencies, we will assume all frequencies to be 1 for the time being (manually added to Excel):
if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Khare2015, here("data_in", "for_func", "Khare2015.csv"))
single_long("Khare2015", "pvdJ spent media", "Not mentioned", "pvdJ spent media", "Ecoli_K12")



# (Tri): Sandberg2017, subsetted by environment. 
## (Tri): Sandberg2017_ac:
Sandberg2017_ac <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_ac_usable.csv"))
Sandberg2017_ac <- Sandberg2017_ac %>% 
  replace(is.na(.), 0) %>%
  filter(details!= "intergenic")
Sandberg2017_ac <- gather(Sandberg2017_ac, population, frequency, "A1 F56 I1 R1" : "A4 F55 I1 R1", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Sandberg2017_ac, here("data_in", "for_func", "Sandberg2017_ac.csv"))

## (Tri): Sandberg2017_glu_ac:
Sandberg2017_glu_ac <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_ac_usable.csv"))
Sandberg2017_glu_ac <- Sandberg2017_glu_ac %>% 
  filter(details!= "intergenic") %>% 
  replace(is.na(.), 0)
### (Tri): transmute() here is used to collapse the clones & replicates from the same population.
Sandberg2017_glu_ac <- Sandberg2017_glu_ac %>% 
  transmute(gene = gene, details = details, 
            A7 = rowSums(Sandberg2017_glu_ac[, 4:15]), A8 = rowSums(Sandberg2017_glu_ac[, 16:28]), 
            A9 = rowSums(Sandberg2017_glu_ac[, 29:38]))
Sandberg2017_glu_ac <- gather(Sandberg2017_glu_ac, population, frequency, "A7" : "A9", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Sandberg2017_glu_ac, here("data_in", "for_func", "Sandberg2017_glu_ac.csv"))

## (Tri): Sandberg2017_glu_gly:
Sandberg2017_glu_gly <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_gly_usable.csv"))
Sandberg2017_glu_gly <- Sandberg2017_glu_gly %>% 
  filter(details!= "intergenic") %>% 
  replace(is.na(.), 0)
Sandberg2017_glu_gly <- Sandberg2017_glu_gly %>% 
  transmute(gene = gene, details = details, 
            A4 = rowSums(Sandberg2017_glu_gly[, 4:7]), A5 = rowSums(Sandberg2017_glu_gly[, 8:12]), 
            A6 = rowSums(Sandberg2017_glu_gly[, 13:16]))
Sandberg2017_glu_gly <- gather(Sandberg2017_glu_gly, population, frequency, "A4" : "A6", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Sandberg2017_glu_gly, here("data_in", "for_func", "Sandberg2017_glu_gly.csv"))

## (Tri): Sandberg2017_glu_xyl:
Sandberg2017_glu_xyl <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_xyl_usable.csv"))
Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl %>% 
  filter(details!= "intergenic") %>% 
  replace(is.na(.), 0)
Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl %>% 
  transmute(gene = gene, details = details, 
            A0 = rowSums(Sandberg2017_glu_xyl[, 3]), A1 = rowSums(Sandberg2017_glu_xyl[, 4:13]), 
            A2 = rowSums(Sandberg2017_glu_xyl[, 14:19]), A3 = rowSums(Sandberg2017_glu_xyl[, 20:ncol(Sandberg2017_glu_xyl)]))
Sandberg2017_glu_xyl <- gather(Sandberg2017_glu_xyl, population, frequency, "A1" : "A3", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Sandberg2017_glu_xyl, here("data_in", "for_func", "Sandberg2017_glu_xyl.csv"))

## (Tri): Sandberg2017_xyl:
Sandberg2017_xyl <- read_csv(here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_xyl_usable.csv"))
Sandberg2017_xyl <- Sandberg2017_xyl %>% 
  filter(details!= "intergenic") %>% 
  replace(is.na(.), 0)
Sandberg2017_xyl <- gather(Sandberg2017_xyl, population, frequency, "A1 F116 I1 R1" : "A4 F113 I1 R1", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Sandberg2017_xyl, here("data_in", "for_func", "Sandberg2017_xyl.csv"))
### (Tri): [Some parallel genes in the "Analysis" file appear multiple times.]

single_long("Sandberg2017_ac", "M9 minimal medium + acetate", "1000", "M9 minimal medium + acetate", "Ecoli_K12")
single_long("Sandberg2017_glu_ac", "M9 minimal medium + glucose + acetate", "650", "M9 minimal medium + glucose + acetate", "Ecoli_K12")
single_long("Sandberg2017_glu_gly", "M9 minimal medium + glucose + glycerol", "1170", "M9 minimal medium + glucose + glycerol", "Ecoli_K12")
single_long("Sandberg2017_glu_xyl", "M9 minimal medium + glucose + xylose", "1180", "M9 minimal medium + glucose + xylose", "Ecoli_K12")
single_long("Sandberg2017_xyl", "M9 minimal medium + xylose", "1000", "M9 minimal medium + xylose", "Ecoli_K12")


# (Tri): Griffith2019, subsetted by pH:
## (Tri): Griffith2019_pH_6.5:
Griffith2019_pH_6.5 <- read_csv(here("data_in", "original & usable", "Griffith2019", "Griffith2019_pH_6.5_usable.csv"))
Griffith2019_pH_6.5 <- Griffith2019_pH_6.5 %>% 
  filter(annotation!= "intergenic") %>% 
  replace(is.na(.), 0)
Griffith2019_pH_6.5 <- gather(Griffith2019_pH_6.5, population, frequency, "C-A1-1" : "C-H5-1", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Griffith2019_pH_6.5, here("data_in", "for_func", "Griffith2019_pH_6.5.csv"))


## (Tri): Griffith2019_pH_8:
Griffith2019_pH_8 <- read_csv(here("data_in", "original & usable", "Griffith2019", "Griffith2019_pH_8_usable.csv"))
Griffith2019_pH_8 <- Griffith2019_pH_8 %>% 
  filter(annotation!= "intergenic") %>% 
  replace(is.na(.), 0)
Griffith2019_pH_8 <- gather(Griffith2019_pH_8, population, frequency, "C-G7-1" : "C-D11-1", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Griffith2019_pH_8, here("data_in", "for_func", "Griffith2019_pH_8.csv"))

single_long("Griffith2019_pH_6.5", "M9 minimal medium + xylose", "1000", "pH 6.5", "Ecoli_K12")
single_long("Griffith2019_pH_8", "M9 minimal medium + xylose", "1000", "pH 8", "Ecoli_K12")



# (Tri): Avrani2017:
Avrani2017 <- read_csv(here("data_in", "original & usable", "Avrani2017", "Avrani2017_usable.csv"))
Avrani2017 <- Avrani2017 %>% 
  filter(details!= "intergenic") %>% 
  replace(is.na(.), 0)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Avrani2017, here("data_in", "for_func", "Avrani2017.csv"))
single_long("Avrani2017", "LB", "11 days", "LB", "Ecoli_K12")



# (Tri): Charusanti2010:
Charusanti2010 <- read_csv(here("data_in", "original & usable", "Charusanti2010", "Charusanti2010_usable.csv"))
Charusanti2010 <- Charusanti2010 %>% 
  filter(details!= "intergenic") %>% 
  replace(is.na(.), 0)
Charusanti2010 <- gather(Charusanti2010, population, frequency, "A1" : "A10", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Charusanti2010, here("data_in", "for_func", "Charusanti2010.csv"))
single_long("Charusanti2010", "M9 minimal medium", "50 days", "M9 minimal medium", "Ecoli_K12")



# (Tri): Ramiro2016:
Ramiro2016 <- read_csv(here("data_in", "original & usable", "Ramiro2016", "Ramiro2016_usable.csv"))
Ramiro2016 <- Ramiro2016 %>% 
  filter(details!= "intergenic") %>% 
  replace(is.na(.), 0)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Ramiro2016, here("data_in", "for_func", "Ramiro2016.csv"))
single_long("Ramiro2016", "RPMI", "60-120", "RPMI", "Ecoli_K12")



# (Tri): Anand2019 (only 1 founder had 3 replicates evolving in the same condition)(gene names after underscores indicate genes missing):
Anand2019_menF_entC_ubiC <- read_csv(here("data_in", "original & usable", "Anand2019", "Anand2019_menF_entC_ubiC_usable.csv"))
Anand2019_menF_entC_ubiC <- Anand2019_menF_entC_ubiC %>% 
  filter(details!= "intergenic") %>% 
  replace(is.na(.), 0)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Anand2019_menF_entC_ubiC, here("data_in", "for_func", "Anand2019_menF_entC_ubiC.csv"))
single_long("Anand2019_menF_entC_ubiC", "M9 minimal medium", "38-47 days", "M9 minimal medium", "Ecoli_K12")



# (Tri): Tenaillon2012:
Tenaillon2012 <- read_csv(here("data_in", "original & usable", "Tenaillon2012", "Tenaillon2012_usable.csv"))
Tenaillon2012 <-Tenaillon2012 %>% 
  filter(Details != "intergenic") %>%
  replace(is.na(.), 0)
Tenaillon2012 <- gather(Tenaillon2012, population, frequency, "A0 F1 I1 R1": "A143 F1 I1 R1", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Tenaillon2012, here("data_in", "for_func", "Tenaillon2012.csv"))
single_long("Tenaillon2012", "Davis minimal medium", "2000", "Davis minimal medium", "Ecoli_K12")
### (Tri): I tried to run Tenaillon2012 with single_wide(), but to no avail.
### (Tri): A likely reason the "undefined columns selected" error keeps popping up when trying to run with single_wide:
### (Tri): The column names have spaces, which makes them syntactically incorrect.
### (Tri): When running the code, something from the inner workings of the functions used to create single_wide() automatically changes column names to syntactically valid names in order to meet the requirements of some deep-level functions.
### (Tri); Then, what happens is later on in single_wide(), the dataframe num_parallel tries to find population columns whose names were entered BEFORE they were turned into syntactically correct names!
### (Tri): Since the population columns now have new, syntax-abiding names, num_parallel can't find the old names, hence the "undefined columns selected" error.
### (Tri): So, to tackle this problem, I use make.names(), which makes syntactically valid names from character vectors:
# colnames(Tenaillon2012) <- make.names(colnames(Tenaillon2012), unique = TRUE)
# if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Tenaillon2012, here("data_in", "for_func", "Tenaillon2012.csv"))
### (Tri): The "undefined columns selected" problem is solved, but the "need finite 'ylim' values" appears at line 101 of single_wide().
# single_wide("Tenaillon2012", colnames(Tenaillon2012)[3:ncol(Tenaillon2012)], "Davis minimal medium", "Ecoli_K12")



# (Working) Wannier2018:
### This paper originally had two clones sequenced for each population. Collapse the two columns and save the resulting dataset. 
Wannier2018 <- read_csv(here("data_in", "original & usable", "Wannier2018", "Wannier2018_usable.csv"))
Wannier2018 <- Wannier2018 %>% 
  filter(Details != "intergenic") %>%
  replace(is.na(.), 0)
Wannier2018 <- Wannier2018 %>%
  transmute(gene = Gene, A1= `A1 F1 I1 R1`+ `A1 F1 I2 R1`, A2 = `A2 F1 I1 R1`+ `A2 F1 I2 R1`, A3 = `A3 F1 I1 R1`+ `A3 F1 I2 R1`, A4 = `A4 F1 I1 R1`+ `A4 F1 I2 R1`, A5 = `A5 F1 I1 R1`+ `A5 F1 I2 R1`, A6 = `A6 F1 I1 R1`+ `A6 F1 I2 R1`, A7 = `A7 F1 I1 R1`+ `A7 F1 I2 R1`, A8 = `A8 F1 I1 R1`+ `A8 F1 I2 R1`, A9 = `A9 F1 I1 R1`+ `A9 F1 I2 R1`, A10 = `A10 F1 I1 R1`+ `A10 F1 I2 R1`,A11 = `A1 F1 I1 R1`+ `A11 F1 I2 R1`, A12 = `A12 F1 I1 R1`+ `A12 F1 I2 R1`, A13 = `A13 F1 I1 R1`+ `A13 F1 I2 R1`, A14 = `A14 F1 I1 R1`+ `A14 F1 I2 R1`) 
### (Tri): MH's code didn't work (see below), so I converted the dataset to single long & analyze w/ single_long().
Wannier2018 <- gather(Wannier2018, population, frequency, "A1" : "A14", factor_key=TRUE)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Wannier2018, here("data_in", "for_func", "Wannier2018.csv"))
single_long("Wannier2018", "Glucose mininal medium", "1100", "Glucose mininal medium", "Ecoli_K12")

# (Tri) MH's code for Wannier2018:
### (Tri) Extract just the population columns from Wannier2018:
# Wannier2018.matrix <- Wannier2018 %>%
#   select(population) %>%
#   as.matrix(.)
### (Tri) Replace all > 1 values with 1:
# Wannier2018.matrix[Wannier2018.matrix > 0] <-1
### (Tri) Recombine the new population columns with genes:
# Wannier2018 <- cbind(gene = Wannier2018$gene, as.data.frame(Wannier2018.matrix))
# if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(Wannier2018, here("data_in", "Wannier2018.csv"))
### (Tri) The single_wide() call isn't working - "Error in plot.window(...) : need finite 'ylim' values". 
### (Tri) Error traces back to line 101 of single_wide.R, "estimate_pa(full_matrix, ndigits = 4, show.plot = T)".
# single_wide("Wannier2018", c("A1" , "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14"), "Glucose Mininal Medium", "Ecoli_K12")



# (Tri) KuzdzalFick2018:
KuzdzalFick2018 <- read_csv(here("data_in", "original & usable", "KuzdzalFick2018", "KuzdzalFick2018_usable.csv"))
KuzdzalFick2018 <- KuzdzalFick2018 %>% 
  filter(annotation != "intergenic") %>%
  replace(is.na(.), 0)

if (length(Sherlock2013_out) > 0) {   Sherlock2013 <- Sherlock2013[-Sherlock2013_out,] } write_csv(KuzdzalFick2018, here("data_in", "for_func", "KuzdzalFick2018.csv"))
single_long("KuzdzalFick2018", "YPD", "28 days", "YPD", "Sac")

####################################

# (Tri) Combine all the analysis files:
file_list <- list.files("data_out/analyses", full.names = TRUE)
### (Tri) Load the library "plyr" for this particular line of code only - calling this library earlier would lead to "unused arguments" error when using here().
### (Tri): I think it's because the function here() is also included in the library "plyr", but it works differently.
master_analyses <- plyr::ldply(file_list, read_csv)
write_csv(master_analyses, "data_out/master_analyses.csv")

####################################
# (Tri) In development: Integrate gather() to single_wide() - easier to process, less repetitive codes.
