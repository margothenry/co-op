#' Calculations for ERConstraint
#'
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for any dataset.
#'
#' @param paper The name of the paper containing the dataset of interest.
#' @param dataset_name The actual name of the dataset (the part before "_usable.csv")
#' @param timepoint_pressure_info Is the dataset for a single timepoint (Value: "single") or multiple timepoints/selective pressures (Value: "multiple")? This information will also be used later to call out different functions according to how the data is organized.
#' @param structure How the data is organized. Values: "wide" & "long". 
###### If the data is for a single timepoint, a wide dataset would contain one column for each population & a long dataset would contain a single column specifying the population.
###### If the data is for multiple timepoints/selective, a wide dataset would contain one column for each timepoint & a long dataset would contain a single column specifying the selective pressure.
#' @param environment The environment in which the experiment occured.
#' @param generations Timepoint(s) in the data, if generations are used to notate. Must be numeric.
#' @param selective_pressure A list of the selective pressures in the data. i.e: temperatures, media, stressors.
#' @param species Specifies if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species when prompted.
#' @param ploidy Haploid, diploid, etc. For E. coli, it's always haploid.
#' (190620: In development) @param collapseMutations Specifies whether to run the analysis at the gene level or on distinct mutations within a gene. The default is at the gene level, i.e. to collapse all different mutations within a gene to one entry in the analysis.
#' @param numgenes The number of genes of the investigated species. If the species specified above is in the database, there's no need to enter a number here.
#' @param strain_info The specifics of the strain (i.e. the "mucoid" in "Wielgoss2016_mucoid")
#' @param population (For single_wide datasets only) A vector containing the names of population columns.
#' @param days Timepoint(s) in the data, if days are used to notate. Remember to call with "days = ". Must be numeric.
#' @param flasks Timepoint(s) in the data, if flasks are are used to notate. Remember to call with "flasks = ". Must be numeric. Only 1 of the 3 potential timepoint types shall be called.
#' @param who_analyzed Who analyzed this dataset? Use initials of 1st & last names.
#' @return A table with all the calculated information.
#' @export
###########################
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
########################### Only the libraries below are needed for the function to work. The ones above are for the data-cleaning steps in the examples.
library(dplyr)
library(sjmisc)
source("R/dgconstraint/functions/single_long.R")
source("R/dgconstraint/functions/single_wide.R")
source("R/dgconstraint/functions/multiple_long.R")
source("R/dgconstraint/functions/multiple_wide.R")
###########################
ERConstraint <- function(paper, dataset_name, timepoint_pressure_info, structure, environment, generations = NA, selective_pressure, species  = NA, ploidy, numgenes = NA, strain_info = NA, population = NA, days = NA, flasks = NA, who_analyzed){

  if(timepoint_pressure_info == "single" & structure == "long"){
    .single_long(paper = paper, dataset_name = dataset_name, environment = environment, generations = generations, selective_pressure = selective_pressure, species = species, ploidy = ploidy, numgenes = numgenes, 
                 strain_info = strain_info, days = days, flasks = flasks, who_analyzed = who_analyzed)
  }
  if (timepoint_pressure_info == "single" & structure == "wide"){
    .single_wide(paper = paper, dataset_name = dataset_name, environment = environment, generations = generations, selective_pressure = selective_pressure, species = species, ploidy = ploidy, numgenes = numgenes, 
                 strain_info = strain_info, population = population, days = days, flasks = flasks, who_analyzed = who_analyzed)
  }
  if (timepoint_pressure_info == "multiple" & structure == "long"){
    .multiple_long(paper = paper, dataset_name = dataset_name, environment = environment, generations = generations, selective_pressure = selective_pressure, species = species, ploidy = ploidy, numgenes = numgenes, 
                   strain_info = strain_info, days = days, flasks = flasks, who_analyzed = who_analyzed)
  }
  if (timepoint_pressure_info == "multiple" & structure == "wide"){
    .multiple_wide(paper = paper, dataset_name = dataset_name, environment = environment, generations = generations, selective_pressure = selective_pressure, species = species, ploidy = ploidy, numgenes = numgenes, 
                   strain_info = strain_info, days = days, flasks = flasks, who_analyzed = who_analyzed)
  }
}

# Examples for the 4 different data types:
## single_long:
Sandra2016 <- read_csv(here("data_in", "original & usable", "Sandra2016", "Sandra2016_usable.csv"))
Sandra2016 <- clean_names(Sandra2016, case = "snake")
colnames(Sandra2016) <- tolower(colnames(Sandra2016))
Sandra2016_out <- c(grep(out_patterns_column_gene, Sandra2016$gene), grep(out_patterns_column_details, Sandra2016$details))
if (length(Sandra2016_out) > 0) {   
  Sandra2016 <- Sandra2016[-Sandra2016_out,] 
} 
Sandra2016 <- Sandra2016 %>%
  select(gene, population, frequency)
Sandra2016 <- Sandra2016 %>%
  replace_na(value = 0)
Sandra2016$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandra2016$gene)

write_csv(Sandra2016, here("data_in", "for_func", "Sandra2016.csv"))
ERConstraint(paper = "Sandra2016", dataset_name = "Sandra2016", timepoint_pressure_info = "single", structure = "long", environment = "M9 minimal medium", days = "8", 
            selective_pressure = "Combined UVA and UVB irradiation", species = "Ecoli_K12", who_analyzed = "TL", ploidy = "haploid")



## single_wide:
### (TL): I honestly could not figure out why single_wide() does not work. The error is at lines 47 & 48 of the script for single_wide ("single_wide.R").
###       Error returns as "Error in UseMethod("drop_na_") : no applicable method for 'drop_na_' applied to an object of class "character"", even though the same drop_na commands are used in the other functions (which also run objects of class "character" through them).



## multiple_long:
Jerison2017 <- read_csv(here("data_in", "original & usable", "Jerison2017", "Jerison2017_usable.csv"))
Jerison2017 <- clean_names(Jerison2017, case = "snake")
colnames(Jerison2017) <- tolower(colnames(Jerison2017))
Jerison2017_out <- c(grep(out_patterns_column_gene, Jerison2017$gene), grep(out_patterns_column_details, Jerison2017$details))
if (length(Jerison2017_out) > 0) {   
  Jerison2017 <- Jerison2017[-Jerison2017_out,] 
} 
Jerison2017 <- Jerison2017 %>%
  select(gene, population, selective_pressure, frequency)
Jerison2017 <- Jerison2017 %>%
  replace_na(value = 0)
Jerison2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Jerison2017$gene)

write_csv(Jerison2017, here("data_in", "for_func", "Jerison2017.csv"))
ERConstraint(paper = "Jerison2017", dataset_name = "Jerison2017", timepoint_pressure_info = "multiple", structure = "long", environment = "SC", generations = "500", selective_pressure = c("HT", "OT"), 
              species = "Sac", who_analyzed = "MH", ploidy = "haploid")



## multiple_wide:
Tenaillon2016 <- read_csv(here("data_in", "original & usable", "Tenaillon2016", "Tenaillon2016_usable.csv"))
Tenaillon2016 <- clean_names(Tenaillon2016, case = "snake")
colnames(Tenaillon2016) <- tolower(colnames(Tenaillon2016))
Tenaillon2016_out <- c(grep(out_patterns_column_gene, Tenaillon2016$gene), grep(out_patterns_column_details, Tenaillon2016$details))
if (length(Tenaillon2016_out) > 0) {
  Tenaillon2016 <- Tenaillon2016[-Tenaillon2016_out,]
}
Tenaillon2016 <- Tenaillon2016 %>% 
  transmute(gene, population, "500" = `x500_i1_r1`+`x500_i2_r1`, "1000" =`x1000_i1_r1`+`x1000_i2_r1`, "1500" =`x1500_i1_r1`+`x1500_i2_r1`,  "2000"= `x2000_i1_r1`+`x2000_i2_r1`, "5000"=`x5000_i1_r1`+`x5000_i2_r1`, "10000"= `x10000_i1_r1`+`x10000_i2_r1`, "15000"=`x15000_i1_r1`+`x15000_i2_r1`,"20000"=`x20000_i1_r1`+`x20000_i2_r1`,"30000"= `x30000_i1_r1`+`x30000_i2_r1`,"40000"=`x40000_i1_r1`+`x40000_i2_r1`,"50000"=`x50000_i1_r1`+`x50000_i2_r1`) 
Tenaillon2016 <- Tenaillon2016 %>%
  replace_na(value = 0)
Tenaillon2016$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Tenaillon2016$gene)

write_csv(Tenaillon2016, here("data_in", "for_func", "Tenaillon2016.csv"))
ERConstraint(paper = "Tenaillon2016", dataset_name = "Tenaillon2016", timepoint_pressure_info = "multiple", structure = "wide", environment = "Davis minimal medium", 
              generations = c("500", "1000", "1500", "2000", '5000', '10000', '15000', '20000', '30000', '40000','50000'), selective_pressure = "Davis minimal medium", 
              species = "Ecoli_K12", who_analyzed = "MH", ploidy = "haploid")

