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
#' @param population_range
#' @param collapseMutations Specifies whether to run the analysis at the gene level or on distinct mutations within a gene. The default is at the gene level, i.e. to collapse all different mutations within a gene to one entry in the analysis.
#' @param numgenes The number of genes of the investigated species. If the species specified above is in the database, there's no need to enter a number here.
#' @return A table with all the calculated information.
#' @export 
#' @examples: single_wide_2
#'
source("R/dgconstraint/functions/single_long.R")

Tenaillon2012 <- read_csv(here("data_in", "original & usable", "Tenaillon2012", "Tenaillon2012_usable.csv"))
Tenaillon2012 <-Tenaillon2012 %>% 
  filter(Details != "intergenic") %>%
  replace(is.na(.), 0)

write_csv(Tenaillon2012, here("data_in", "for_func", "Tenaillon2012.csv"))

single_wide_2 <- function(paper, environment, generation, selective_pressure, species = NA, population_range, collapseMutations = TRUE, numgenes = NA){
  geneNumbers <- read_csv(file.path(getwd(),"R/dgconstraint/inst/GeneDatabase.csv"), col_types = cols())
  
  data <- read_csv(file.path(getwd(), "data_in", "for_func", paste0(paper, ".csv")), col_types = cols())
  
  # (Tri) If the species name is found in the database, the "numgenes" can be retrieved there.
  if (species %in% geneNumbers$Species){
    numgenes <- filter(geneNumbers, Species == species)$NumGenes  
  }
  
  # (Tri) Otherwise, enter "numgenes" by hand.
  if(is.na(numgenes)){
    prompt <- "Your species is unspecified or not in our database. How many genes does it have? \n"
    numgenes <- as.numeric(readline(prompt))
  }
  # (Tri) Convert column range to numeric form:
  cols.num <- population_range
  paper[cols.num,1] <- sapply(paper[cols.num,1], as.numeric)
  cols.num <- paper[cols.num,1]
  # (Tri) gather() changes data from wide to long format:
  paper <- gather(paper, population, frequency, population_range, factor_key=TRUE)
  single_long(paper, environment, generation, selective_pressure, species = NA, numgenes = NA)
}

Tenaillon2012 <- read_csv(here("data_in", "original & usable", "Tenaillon2012", "Tenaillon2012_usable.csv"))
Tenaillon2012 <-Tenaillon2012 %>% 
  filter(Details != "intergenic") %>%
  replace(is.na(.), 0)

write_csv(Tenaillon2012, here("data_in", "for_func", "Tenaillon2012.csv"))
single_wide_2(Tenaillon2012, "Davis minimal medium", "2000", "Davis minimal medium", "Ecoli_K12", c(3:115))
  
  