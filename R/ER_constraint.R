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
if( !require( "pacman" ) ) { install.packages( "pacman" ) }
pacman::p_load(
  here,
  tidyverse,
  here,
  dgconstraint,
  R.utils,
  janitor,
  readr,
  devtools,
  Hmisc,
  naniar
)
########################### Only the libraries below are needed for the function to work. The ones above are for the data-cleaning steps in the examples.

if( !require( "pacman" ) ) { install.packages( "pacman" ) }
pacman::p_load(
  sjmisc
)
source("R/dgconstraint/functions/single_long.R")
source("R/dgconstraint/functions/multiple_long.R")
source("R/dgconstraint/functions/multiple_wide.R")

###########################

ERConstraint <- function(
  paper,
  dataset_name,
  timepoint_pressure_info,
  structure,
  environment,
  generations = NA,
  selective_pressure,
  species = NA,
  ploidy,
  numgenes = NA,
  strain_info = NA,
  population = NA,
  days = NA,
  flasks = NA,
  who_analyzed
  )
{
  
  if(timepoint_pressure_info == "single" & structure == "long"){
    .single_long(
      paper = paper,
      dataset_name = dataset_name,
      environment = environment,
      generations = generations,
      selective_pressure = selective_pressure,
      species = species,
      ploidy = ploidy,
      numgenes = numgenes, 
      strain_info = strain_info,
      days = days,
      flasks = flasks,
      who_analyzed = who_analyzed
      )
  }
  if (timepoint_pressure_info == "multiple" & structure == "long"){
    .multiple_long(
      paper = paper,
      dataset_name = dataset_name,
      environment = environment,
      generations = generations,
      selective_pressure = selective_pressure,
      species = species,
      ploidy = ploidy,
      numgenes = numgenes, 
      strain_info = strain_info,
      days = days,
      flasks = flasks, 
      who_analyzed = who_analyzed
      )
  }
  if (timepoint_pressure_info == "multiple" & structure == "wide"){
    .multiple_wide(
      paper = paper,
      dataset_name = dataset_name,
      environment = environment,
      generations = generations,
      selective_pressure = selective_pressure,
      species = species,
      ploidy = ploidy, 
      numgenes = numgenes, 
      strain_info = strain_info,
      days = days, 
      flasks = flasks, 
      who_analyzed = who_analyzed
      )
  }
}
