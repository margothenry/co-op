#' IN DEVELOPMENT
#' The difference between this & single_wide() is this aims to gather the wide data into long format & use single_long() to run the analysis.
#' Calculations for a Single Wide Dataset
#'
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a dataset with a single generation.
#' 
#' @param paper The data in .csv that you want to analyze.
#' @param environment The environment in which the experiment occured.
#' @param generation The generation the sequencing took place. Could also be a timepoint if the generation isn't specified in the paper. Make sure to include "units" (e.g. days, flasks) for non-generation entries.
#' @param selective_pressure A list of the selective pressures in the data. i.e: temperatures, media, stressors.
#' @param species Specifies if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species when prompted.
#' @param collapseMutations Specifies whether to run the analysis at the gene level or on distinct mutations within a gene. The default is at the gene level, i.e. to collapse all different mutations within a gene to one entry in the analysis.
#' @param numgenes The number of genes of the investigated species. If the species specified above is in the database, there's no need to enter a number here.
#' @return A table with all the calculated information.
#' @export 
#' @examples: single_wide_2
#'
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
source("R/dgconstraint/functions/single_long.R")

single_wide_2 <- function(paper, dataset_name, environment, generations, days, flasks, selective_pressure, species = NA, ploidy, numgenes = NA, strain_info){
  paper <- gather(paper, population, frequency, 3:ncol(paper), factor_key=TRUE)
  single_long(paper, dataset_name, environment, generations, days, flasks, selective_pressure, species, ploidy, numgenes, strain_info)
}

Tenaillon2012 <- read_csv(here("data_in", "original & usable", "Tenaillon2012", "Tenaillon2012_usable.csv"))
Tenaillon2012 <- clean_names(Tenaillon2012, case = "snake")
colnames(Tenaillon2012) <- tolower(colnames(Tenaillon2012))
Tenaillon2012 <- Tenaillon2012 %>%
  replace(is.na(.), 0) %>% 
  filter(details != "intergenic") 
Tenaillon2012_out <- c(grep(",", Tenaillon2012$gene), grep("/", Tenaillon2012$gene), grep("[*]", Tenaillon2012$gene),  
                       grep("-", Tenaillon2012$gene))
if (length(Tenaillon2012_out) > 0) {   
  Tenaillon2012 <- Tenaillon2012[-Tenaillon2012_out,] 
} 

write_csv(Tenaillon2012, here("data_in", "for_func", "Tenaillon2012.csv"))
single_wide_2(paper = "Tenaillon2012", dataset_name = "Tenaillon2012", environment = "Davis minimal medium", generations = "2000", 
            selective_pressure = "heat (42.2 degrees C)", species = "Ecoli_K12", ploidy = "haploid")  
### (190712) No applicable method for 'gather_' applied to an object of class "character"?
### (190712) Call out col name by func?
  