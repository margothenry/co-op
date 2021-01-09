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

out_patterns_column_gene <- c(",|/|\\[|-|\\b0\\b")
out_patterns_column_details <- c("prophage|extragenic|upstream|intergenic|5' UTR|LTR|Intron")

# Anand2019
Anand2019_menF_entC_ubiC <- read_csv(
  here("data_in", "original & usable", "Anand2019", "Anand2019_menF_entC_ubiC_usable.csv")
  )

Anand2019_menF_entC_ubiC <- clean_names(Anand2019_menF_entC_ubiC, case = "snake")

colnames(Anand2019_menF_entC_ubiC) <- tolower(colnames(Anand2019_menF_entC_ubiC))

Anand2019_menF_entC_ubiC_out <- c(
  grep(out_patterns_column_gene, Anand2019_menF_entC_ubiC$gene),
  grep(out_patterns_column_details, Anand2019_menF_entC_ubiC$details)
  )

if (length(Anand2019_menF_entC_ubiC_out) > 0) {   
  Anand2019_menF_entC_ubiC <- Anand2019_menF_entC_ubiC[-Anand2019_menF_entC_ubiC_out,] 
} 

Anand2019_menF_entC_ubiC <- Anand2019_menF_entC_ubiC %>%
  select(gene, population, frequency)

Anand2019_menF_entC_ubiC <- Anand2019_menF_entC_ubiC %>%
  replace_na(value = 0)

Anand2019_menF_entC_ubiC$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Anand2019_menF_entC_ubiC$gene)

write_csv(Anand2019_menF_entC_ubiC, here("data_in", "for_func", "Anand2019_menF_entC_ubiC.csv"))

ERConstraint(
  paper = "Anand2019",
  dataset_name = "Anand2019_menF_entC_ubiC",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  days = "41", 
  selective_pressure = "iron (FeSO4)",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Missing menF, entC, & ubiC genes"
  ) 

# Avrani2017
Avrani2017 <- read_csv(
  here("data_in", "original & usable", "Avrani2017", "Avrani2017_usable.csv")
  )

Avrani2017 <- clean_names(Avrani2017, case = "snake")

colnames(Avrani2017) <- tolower(colnames(Avrani2017))

Avrani2017_out <- c(
  grep(out_patterns_column_gene, Avrani2017$gene),
  grep(out_patterns_column_details, Avrani2017$details)
  )

if (length(Avrani2017_out) > 0) {   
  Avrani2017 <- Avrani2017[-Avrani2017_out,] 
} 

Avrani2017 <- Avrani2017 %>%
  select(gene, population, frequency)

Avrani2017 <- Avrani2017 %>%
  replace_na(value = 0)

Avrani2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Avrani2017$gene)

write_csv(Avrani2017, here("data_in", "for_func", "Avrani2017.csv"))

ERConstraint(
  paper = "Avrani2017",
  dataset_name = "Avrani2017",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "LB",
  days = "11", 
  selective_pressure = "Resource exhaustion",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Boyle2017
Boyle2017 <- read_csv(
  here("data_in", "original & usable", "Boyle2017", "Boyle2017_usable.csv")
  )

Boyle2017 <- clean_names(Boyle2017, case = "snake")

colnames(Boyle2017) <- tolower(colnames(Boyle2017))

Boyle2017_out <- c(
  grep(out_patterns_column_gene, Boyle2017$gene),
  grep(out_patterns_column_details, Boyle2017$details)
  )

if (length(Boyle2017_out) > 0) {   
  Boyle2017 <- Boyle2017[-Boyle2017_out,] 
} 

Boyle2017 <- Boyle2017 %>%
  select(gene, population, frequency)

Boyle2017 <- Boyle2017 %>%
  replace_na(value = 0)

Boyle2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Boyle2017$gene)

write_csv(Boyle2017, here("data_in", "for_func", "Boyle2017.csv"))

ERConstraint(
  paper = "Boyle2017",
  dataset_name = "Boyle2017",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "LB",
  generations = "49", 
  selective_pressure = "LB",
  species = "P_aeruginosa_PA14",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "cbrA deleted"
  )

# Chib2017
Chib2017 <- read_csv(
  here("data_in", "original & usable", "Chib2017", "Chib2017_usable.csv")
  )

Chib2017 <- clean_names(Chib2017, case = "snake")

colnames(Chib2017) <- tolower(colnames(Chib2017))

Chib2017_out <- c(
  grep(out_patterns_column_gene, Chib2017$gene),
  grep(out_patterns_column_details, Chib2017$details)
  )

if (length(Chib2017_out) > 0) {   
  Chib2017 <- Chib2017[-Chib2017_out,] 
}

Chib2017 <- Chib2017 %>%
  select(gene, population, frequency)

Chib2017 <- Chib2017 %>%
  replace_na(value = 0)

Chib2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Chib2017$gene)

write_csv(Chib2017, here("data_in", "for_func", "Chib2017.csv"))

ERConstraint(
  paper = "Chib2017",
  dataset_name = "Chib2017",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "LB",
  days = "28", 
  selective_pressure = "Prolonged stationary phase",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Conrad2009
## Conrad2009_A_to_E:
Conrad2009_A_to_E <- read_csv(
  here("data_in", "original & usable", "Conrad2009", "Conrad2009_A_to_E_usable.csv")
  )

Conrad2009_A_to_E <- clean_names(Conrad2009_A_to_E, case = "snake")

colnames(Conrad2009_A_to_E) <- tolower(colnames(Conrad2009_A_to_E))

Conrad2009_A_to_E_out <- c(
  grep(out_patterns_column_gene, Conrad2009_A_to_E$gene),
  grep(out_patterns_column_details, Conrad2009_A_to_E$details)
  )

if (length(Conrad2009_A_to_E_out) > 0) {   
  Conrad2009_A_to_E <- Conrad2009_A_to_E[-Conrad2009_A_to_E_out,] 
} 

Conrad2009_A_to_E <- Conrad2009_A_to_E %>%
  select(gene, population, frequency)

Conrad2009_A_to_E <- Conrad2009_A_to_E %>%
  replace_na(value = 0)

Conrad2009_A_to_E$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Conrad2009_A_to_E$gene)

write_csv(Conrad2009_A_to_E, here("data_in", "for_func", "Conrad2009_A_to_E.csv"))

ERConstraint(
  paper = "Conrad2009",
  dataset_name = "Conrad2009_A_to_E",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9",
  generations = "1100", 
  selective_pressure = "Lactate",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

## Conrad2009_F_to_K:
Conrad2009_F_to_K <- read_csv(
  here("data_in", "original & usable", "Conrad2009", "Conrad2009_F_to_K_usable.csv")
  )

Conrad2009_F_to_K <- clean_names(Conrad2009_F_to_K, case = "snake")

colnames(Conrad2009_F_to_K) <- tolower(colnames(Conrad2009_F_to_K))

Conrad2009_F_to_K_out <- c(
  grep(out_patterns_column_gene, Conrad2009_F_to_K$gene),
  grep(out_patterns_column_details, Conrad2009_F_to_K$details)
  )

if (length(Conrad2009_F_to_K_out) > 0) {   
  Conrad2009_F_to_K <- Conrad2009_F_to_K[-Conrad2009_F_to_K_out,] 
} 

Conrad2009_F_to_K <- Conrad2009_F_to_K %>%
  select(gene, population, frequency)

Conrad2009_F_to_K <- Conrad2009_F_to_K %>%
  replace_na(value = 0)

Conrad2009_F_to_K$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Conrad2009_F_to_K$gene)

write_csv(Conrad2009_F_to_K, here("data_in", "for_func", "Conrad2009_F_to_K.csv"))

ERConstraint(
  paper = "Conrad2009",
  dataset_name = "Conrad2009_F_to_K",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9",
  generations = "750", 
  selective_pressure = "Lactate",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Creamer2016 
Creamer2016 <- read_csv(
  here("data_in", "original & usable", "Creamer2016", "Creamer2016_usable.csv")
  )

Creamer2016 <- clean_names(Creamer2016, case = "snake")

colnames(Creamer2016) <- tolower(colnames(Creamer2016))

Creamer2016_out <- c(
  grep(out_patterns_column_gene, Creamer2016$gene),
  grep(out_patterns_column_details, Creamer2016$details)
  )

if (length(Creamer2016_out) > 0) {
  Creamer2016 <- Creamer2016[-Creamer2016_out,] 
} 

Creamer2016 <- Creamer2016 %>% 
  transmute(
    gene = gene, details = details, 
    "a1" =  a1_1,
    'a3' = a3_1,
    "b1" = b1_1,
    "c1" = c1_1,
    c3 = rowSums(Creamer2016[, 5:6], na.rm = TRUE),
    "d5" = d5_1,
    a5 = rowSums(Creamer2016[, 8:9], na.rm = TRUE), 
    e1 = rowSums(Creamer2016[, 10:11], na.rm = TRUE),
    "h1" = h1_1,
    "h3" = h3_1,
    "g3" = g3_1,
    g5 = rowSums(Creamer2016[, 15:16], na.rm = TRUE)
  )

Creamer2016 <- gather(
  Creamer2016,
  population,
  frequency,
  "a1" : "g5",
  factor_key = TRUE )

Creamer2016 <- Creamer2016 %>%
  select(gene, population, frequency)

Creamer2016 <- Creamer2016 %>%
  replace_na(value = 0)

Creamer2016$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Creamer2016$gene)

write_csv(Creamer2016, here("data_in", "for_func", "Creamer2016.csv"))

ERConstraint(
  paper = "Creamer2016",
  dataset_name = "Creamer2016",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "LBK medium",
  generations = "2000", 
  selective_pressure = "benzoate",
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid"
  )

# Dai2018
Dai2018 <- read_csv(
  here("data_in", "original & usable", "Dai2018", "Dai2018_usable.csv")
  )

Dai2018 <- clean_names(Dai2018, case = "snake")

colnames(Dai2018) <- tolower(colnames(Dai2018))

Dai2018_out <- c(
  grep(out_patterns_column_gene, Dai2018$gene),
  grep(out_patterns_column_details, Dai2018$details)
  )

if (length(Dai2018_out) > 0) {   
  Dai2018 <- Dai2018[-Dai2018_out,] 
}

Dai2018 <- Dai2018 %>%
  select(gene, population, frequency)

Dai2018 <- Dai2018 %>%
  replace_na(value = 0)

Dai2018$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Dai2018$gene)

write_csv(Dai2018, here("data_in", "for_func", "Dai2018.csv"))

ERConstraint(
  paper = "Dai2018",
  dataset_name = "Dai2018",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Minimal glucose medium",
  days = "40", 
  selective_pressure = "Minimal glucose medium",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Deatherage2017
Deatherage2017 <- read_csv(
  here("data_in", "original & usable", "Deatherage2017", "Deatherage2017_usable.csv")
  )

Deatherage2017 <- clean_names(Deatherage2017, case = "snake")

colnames(Deatherage2017) <- tolower(colnames(Deatherage2017))

Deatherage2017_out <- c(
  grep(out_patterns_column_gene, Deatherage2017$gene),
  grep(out_patterns_column_details, Deatherage2017$details)
  )

if (length(Deatherage2017_out) > 0) {   
  Deatherage2017 <- Deatherage2017[-Deatherage2017_out,] 
} 

Deatherage2017 = Deatherage2017 %>% tidyr::replace_na(
  list(
    f1_i1_r1 = 0,
    f1_i2_r1 = 0,
    f1_i3_r1 = 0
  )
) 

Deatherage2017 <- Deatherage2017 %>%
  transmute(
    gene = gene, 
    population = population, 
    "frequency" = `f1_i1_r1` + `f1_i2_r1` + `f1_i3_r1`
    )

Deatherage2017 <- Deatherage2017 %>%
  replace_na(value = 0)

Deatherage2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Deatherage2017$gene)

write_csv(Deatherage2017, here("data_in", "for_func", "Deatherage2017.csv"))

ERConstraint(
  paper = "Deatherage2017",
  dataset_name = "Deatherage2017",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Glucose minimal medium",
  generations = "2000", 
  selective_pressure = "Temperature",
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid"
  )

# Du2019
Du2019 <- read_csv(
  here("data_in", "original & usable", "Du2019", "Du2019_usable.csv")
  )

Du2019 <- clean_names(Du2019, case = "snake")

colnames(Du2019) <- tolower(colnames(Du2019))

names(Du2019) <- gsub("x", "", names(Du2019))

Du2019_out <- c(
  grep(out_patterns_column_gene, Du2019$gene),
  grep(out_patterns_column_details, Du2019$details)
  )

if (length(Du2019_out) > 0) {   
  Du2019 <- Du2019[-Du2019_out,] 
} 

Du2019 <- spread(Du2019, flask, value)

Du2019 <- gather(Du2019, population, frequency, "aa1" : "aa6", factor_key=TRUE)

Du2019 <- Du2019 %>%
  select(gene, population, "intermediate", "late")

Du2019 <- Du2019 %>%
  replace_na(value = 0)

Du2019$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Du2019$gene)

write_csv(Du2019, here("data_in", "for_func", "Du2019.csv"))

ERConstraint(
  paper = "Du2019",
  dataset_name = "Du2019",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "Minimal glucose medium",
  flasks = c("intermediate", "late"), 
  selective_pressure = "pH 5.5",
  species = "Ecoli_K12",
  who_analyzed = "L",
  ploidy = "haploid"
  )

# Echenique2019
## Echenique2019_ADE2:
Echenique2019_ADE2 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_ADE2_usable.csv")
  )

Echenique2019_ADE2 <- clean_names(Echenique2019_ADE2, case = "snake")

colnames(Echenique2019_ADE2) <- tolower(colnames(Echenique2019_ADE2))

Echenique2019_ADE2_out <- c(
  grep(out_patterns_column_gene, Echenique2019_ADE2$gene),
  grep(out_patterns_column_details, Echenique2019_ADE2$details)
  )

if (length(Echenique2019_ADE2_out) > 0) {   
  Echenique2019_ADE2 <- Echenique2019_ADE2[-Echenique2019_ADE2_out,] 
} 

Echenique2019_ADE2 <- Echenique2019_ADE2 %>%
  select(gene, population, frequency)

Echenique2019_ADE2 <- Echenique2019_ADE2 %>%
  replace_na(value = 0)

Echenique2019_ADE2$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_ADE2$gene)

write_csv(Echenique2019_ADE2, here("data_in", "for_func", "Echenique2019_ADE2.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_ADE2",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "ADE2 deleted"
  )

# Echenique2019_BMH1:
Echenique2019_BMH1 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_BMH1_usable.csv")
  )

Echenique2019_BMH1 <- clean_names(Echenique2019_BMH1, case = "snake")

colnames(Echenique2019_BMH1) <- tolower(colnames(Echenique2019_BMH1))

Echenique2019_BMH1_out <- c(
  grep(out_patterns_column_gene, Echenique2019_BMH1$gene),
  grep(out_patterns_column_details, Echenique2019_BMH1$details)
  )

if (length(Echenique2019_BMH1_out) > 0) {   
  Echenique2019_BMH1 <- Echenique2019_BMH1[-Echenique2019_BMH1_out,] 
} 

Echenique2019_BMH1 <- Echenique2019_BMH1 %>%
  select(gene, population, frequency)

Echenique2019_BMH1 <- Echenique2019_BMH1 %>%
  replace_na(value = 0)

Echenique2019_BMH1$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_BMH1$gene)

write_csv(Echenique2019_BMH1, here("data_in", "for_func", "Echenique2019_BMH1.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_BMH1",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "BMH1 deleted"
  )


## Echenique2019_COQ2:
Echenique2019_COQ2 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_COQ2_usable.csv")
  )

Echenique2019_COQ2 <- clean_names(Echenique2019_COQ2, case = "snake")

colnames(Echenique2019_COQ2) <- tolower(colnames(Echenique2019_COQ2))

Echenique2019_COQ2_out <- c(
  grep(out_patterns_column_gene, Echenique2019_COQ2$gene),
  grep(out_patterns_column_details, Echenique2019_COQ2$details)
  )

if (length(Echenique2019_COQ2_out) > 0) {   
  Echenique2019_COQ2 <- Echenique2019_COQ2[-Echenique2019_COQ2_out,] 
} 

Echenique2019_COQ2 <- Echenique2019_COQ2 %>%
  select(gene, population, frequency)

Echenique2019_COQ2 <- Echenique2019_COQ2 %>%
  replace_na(value = 0)

Echenique2019_COQ2$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_COQ2$gene)

write_csv(Echenique2019_COQ2, here("data_in", "for_func", "Echenique2019_COQ2.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_COQ2",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "COQ2 deleted"
  )


## Echenique2019_COX6:
Echenique2019_COX6 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_COX6_usable.csv")
  )

Echenique2019_COX6 <- clean_names(Echenique2019_COX6, case = "snake")

colnames(Echenique2019_COX6) <- tolower(colnames(Echenique2019_COX6))

Echenique2019_COX6_out <- c(
  grep(out_patterns_column_gene, Echenique2019_COX6$gene),
  grep(out_patterns_column_details, Echenique2019_COX6$details)
  )

if (length(Echenique2019_COX6_out) > 0) {   
  Echenique2019_COX6 <- Echenique2019_COX6[-Echenique2019_COX6_out,] 
} 

Echenique2019_COX6 <- Echenique2019_COX6 %>%
  select(gene, population, frequency)

Echenique2019_COX6 <- Echenique2019_COX6 %>%
  replace_na(value = 0)

Echenique2019_COX6$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_COX6$gene)

write_csv(Echenique2019_COX6, here("data_in", "for_func", "Echenique2019_COX6.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_COX6",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "COX6 deleted"
  )


## Echenique2019_CTF19:
Echenique2019_CTF19 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_CTF19_usable.csv")
  )

Echenique2019_CTF19 <- clean_names(Echenique2019_CTF19, case = "snake")

colnames(Echenique2019_CTF19) <- tolower(colnames(Echenique2019_CTF19))

Echenique2019_CTF19_out <- c(
  grep(out_patterns_column_gene, Echenique2019_CTF19$gene),
  grep(out_patterns_column_details, Echenique2019_CTF19$details)
  )

if (length(Echenique2019_CTF19_out) > 0) {   
  Echenique2019_CTF19 <- Echenique2019_CTF19[-Echenique2019_CTF19_out,] 
} 

Echenique2019_CTF19 <- Echenique2019_CTF19 %>%
  select(gene, population, frequency)

Echenique2019_CTF19 <- Echenique2019_CTF19 %>%
  replace_na(value = 0)

Echenique2019_CTF19$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_CTF19$gene)

write_csv(Echenique2019_CTF19, here("data_in", "for_func", "Echenique2019_CTF19.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_CTF19",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "CTF19 deleted"
  )

## Echenique2019_ELP4:
Echenique2019_ELP4 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_ELP4_usable.csv")
  )

Echenique2019_ELP4 <- clean_names(Echenique2019_ELP4, case = "snake")

colnames(Echenique2019_ELP4) <- tolower(colnames(Echenique2019_ELP4))

Echenique2019_ELP4_out <- c(
  grep(out_patterns_column_gene, Echenique2019_ELP4$gene),
  grep(out_patterns_column_details, Echenique2019_ELP4$details)
  )

if (length(Echenique2019_ELP4_out) > 0) {   
  Echenique2019_ELP4 <- Echenique2019_ELP4[-Echenique2019_ELP4_out,] 
} 

Echenique2019_ELP4 <- Echenique2019_ELP4 %>%
  select(gene, population, frequency)

Echenique2019_ELP4 <- Echenique2019_ELP4 %>%
  replace_na(value = 0)

Echenique2019_ELP4$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_ELP4$gene)

write_csv(Echenique2019_ELP4, here("data_in", "for_func", "Echenique2019_ELP4.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_ELP4",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "ELP4 deleted"
  )


## Echenique2019_HXK2:
Echenique2019_HXK2 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_HXK2_usable.csv")
  )

Echenique2019_HXK2 <- clean_names(Echenique2019_HXK2, case = "snake")

colnames(Echenique2019_HXK2) <- tolower(colnames(Echenique2019_HXK2))

Echenique2019_HXK2_out <- c(
  grep(out_patterns_column_gene, Echenique2019_HXK2$gene),
  grep(out_patterns_column_details, Echenique2019_HXK2$details)
  )

if (length(Echenique2019_HXK2_out) > 0) {   
  Echenique2019_HXK2 <- Echenique2019_HXK2[-Echenique2019_HXK2_out,] 
} 

Echenique2019_HXK2 <- Echenique2019_HXK2 %>%
  select(gene, population, frequency)

Echenique2019_HXK2 <- Echenique2019_HXK2 %>%
  replace_na(value = 0)

Echenique2019_HXK2$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_HXK2$gene)

write_csv(Echenique2019_HXK2, here("data_in", "for_func", "Echenique2019_HXK2.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_HXK2",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "HXK2 deleted"
  )


## Echenique2019_IML3:
Echenique2019_IML3 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_IML3_usable.csv")
  )

Echenique2019_IML3 <- clean_names(Echenique2019_IML3, case = "snake")

colnames(Echenique2019_IML3) <- tolower(colnames(Echenique2019_IML3))

Echenique2019_IML3_out <- c(
  grep(out_patterns_column_gene, Echenique2019_IML3$gene),
  grep(out_patterns_column_details, Echenique2019_IML3$details)
  )

if (length(Echenique2019_IML3_out) > 0) {   
  Echenique2019_IML3 <- Echenique2019_IML3[-Echenique2019_IML3_out,] 
} 

Echenique2019_IML3 <- Echenique2019_IML3 %>%
  select(gene, population, frequency)

Echenique2019_IML3 <- Echenique2019_IML3 %>%
  replace_na(value = 0)

Echenique2019_IML3$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_IML3$gene)

write_csv(Echenique2019_IML3, here("data_in", "for_func", "Echenique2019_IML3.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_IML3",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "IML3 deleted"
  )


## Echenique2019_KTI12:
Echenique2019_KTI12 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_KTI12_usable.csv")
  )

Echenique2019_KTI12 <- clean_names(Echenique2019_KTI12, case = "snake")

colnames(Echenique2019_KTI12) <- tolower(colnames(Echenique2019_KTI12))

Echenique2019_KTI12_out <- c(
  grep(out_patterns_column_gene, Echenique2019_KTI12$gene),
  grep(out_patterns_column_details, Echenique2019_KTI12$details)
  )

if (length(Echenique2019_KTI12_out) > 0) {   
  Echenique2019_KTI12 <- Echenique2019_KTI12[-Echenique2019_KTI12_out,] 
} 

Echenique2019_KTI12 <- Echenique2019_KTI12 %>%
  select(gene, population, frequency)

Echenique2019_KTI12 <- Echenique2019_KTI12 %>%
  replace_na(value = 0)

Echenique2019_KTI12$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_KTI12$gene)

write_csv(Echenique2019_KTI12, here("data_in", "for_func", "Echenique2019_KTI12.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_KTI12",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "KTI12 deleted"
  )


## Echenique2019_SOK2:
Echenique2019_SOK2 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_SOK2_usable.csv")
  )

Echenique2019_SOK2 <- clean_names(Echenique2019_SOK2, case = "snake")

colnames(Echenique2019_SOK2) <- tolower(colnames(Echenique2019_SOK2))

Echenique2019_SOK2_out <- c(
  grep(out_patterns_column_gene, Echenique2019_SOK2$gene),
  grep(out_patterns_column_details, Echenique2019_SOK2$details)
  )

if (length(Echenique2019_SOK2_out) > 0) {   
  Echenique2019_SOK2 <- Echenique2019_SOK2[-Echenique2019_SOK2_out,] 
} 

Echenique2019_SOK2 <- Echenique2019_SOK2 %>%
  select(gene, population, frequency)

Echenique2019_SOK2 <- Echenique2019_SOK2 %>%
  replace_na(value = 0)

Echenique2019_SOK2$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_SOK2$gene)

write_csv(Echenique2019_SOK2, here("data_in", "for_func", "Echenique2019_SOK2.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_SOK2",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "SOK2 deleted")


## Echenique2019_VPS29:
Echenique2019_VPS29 <- read_csv(
  here("data_in", "original & usable", "Echenique2019", "Echenique2019_VPS29_usable.csv")
  )

Echenique2019_VPS29 <- clean_names(Echenique2019_VPS29, case = "snake")

colnames(Echenique2019_VPS29) <- tolower(colnames(Echenique2019_VPS29))

Echenique2019_VPS29_out <- c(
  grep(out_patterns_column_gene, Echenique2019_VPS29$gene),
  grep(out_patterns_column_details, Echenique2019_VPS29$details)
  )

if (length(Echenique2019_VPS29_out) > 0) {   
  Echenique2019_VPS29 <- Echenique2019_VPS29[-Echenique2019_VPS29_out,] 
}

Echenique2019_VPS29 <- Echenique2019_VPS29 %>%
  select(gene, population, frequency)

Echenique2019_VPS29 <- Echenique2019_VPS29 %>%
  replace_na(value = 0)

Echenique2019_VPS29$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Echenique2019_VPS29$gene)


write_csv(Echenique2019_VPS29, here("data_in", "for_func", "Echenique2019_VPS29.csv"))

ERConstraint(
  paper = "Echenique2019",
  dataset_name = "Echenique2019_VPS29",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "VPS29 deleted")

# Flynn2014
## Flynn2014_biofilm:
Flynn2014_biofilm <- read_csv(
  here("data_in", "original & usable", "Flynn2014", "Flynn2014_biofilm_usable.csv")
  )

Flynn2014_biofilm <- clean_names(Flynn2014_biofilm, case = "snake")

colnames(Flynn2014_biofilm) <- tolower(colnames(Flynn2014_biofilm))

names(Flynn2014_biofilm) <- gsub("x", "", names(Flynn2014_biofilm))

Flynn2014_biofilm_out <- c(
  grep(out_patterns_column_gene, Flynn2014_biofilm$gene),
  grep(out_patterns_column_details, Flynn2014_biofilm$details)
  )

if (length(Flynn2014_biofilm_out) > 0) {   
  Flynn2014_biofilm <- Flynn2014_biofilm[-Flynn2014_biofilm_out,] 
} 

Flynn2014_biofilm <- Flynn2014_biofilm %>%
  select(gene, population, "102", "150", "264", "396", "450", "540")

Flynn2014_biofilm <- Flynn2014_biofilm %>%
  replace_na(value = 0)

Flynn2014_biofilm$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Flynn2014_biofilm$gene)

write_csv(Flynn2014_biofilm, here("data_in", "for_func", "Flynn2014_biofilm.csv"))

ERConstraint(
  paper = "Flynn2014",
  dataset_name = "Flynn2014_biofilm",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "M63",
  generations = c("102", "150", "264", "396", "450", "540"), 
  selective_pressure = "Polystyrene beads",
  species = "P_aeruginosa_PA14",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "biofilm-evolved"
  )

## Flynn2014_planktonic:
Flynn2014_planktonic <- read_csv(
  here("data_in", "original & usable", "Flynn2014", "Flynn2014_planktonic_usable.csv")
  )

Flynn2014_planktonic <- clean_names(Flynn2014_planktonic, case = "snake")

colnames(Flynn2014_planktonic) <- tolower(colnames(Flynn2014_planktonic))

names(Flynn2014_planktonic) <- gsub("x", "", names(Flynn2014_planktonic))

Flynn2014_planktonic_out <- c(
  grep(out_patterns_column_gene, Flynn2014_planktonic$gene),
  grep(out_patterns_column_details, Flynn2014_planktonic$details)
  )

if (length(Flynn2014_planktonic_out) > 0) {   
  Flynn2014_planktonic <- Flynn2014_planktonic[-Flynn2014_planktonic_out,] 
} 

Flynn2014_planktonic <- Flynn2014_planktonic %>%
  select(gene, population, "396", "540")

Flynn2014_planktonic <- Flynn2014_planktonic %>%
  replace_na(value = 0)

Flynn2014_planktonic$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Flynn2014_planktonic$gene)

write_csv(Flynn2014_planktonic, here("data_in", "for_func", "Flynn2014_planktonic.csv"))

ERConstraint(
  paper = "Flynn2014",
  dataset_name = "Flynn2014_planktonic",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "M63",
  generations = c("396", "540"), 
  selective_pressure = "No polystyrene beads",
  species = "P_aeruginosa_PA14",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "planktonic-evolved"
  )

# Griffith2019   
## Griffith2019_pH_6.5:
Griffith2019_pH_6.5 <- read_csv(
  here("data_in", "original & usable", "Griffith2019", "Griffith2019_pH_6.5_usable.csv")
  )

Griffith2019_pH_6.5 <- clean_names(Griffith2019_pH_6.5, case = "snake")

colnames(Griffith2019_pH_6.5) <- tolower(colnames(Griffith2019_pH_6.5))

Griffith2019_pH_6.5_out <- c(
  grep(out_patterns_column_gene, Griffith2019_pH_6.5$gene),
  grep(out_patterns_column_details, Griffith2019_pH_6.5$details)
  )

if (length(Griffith2019_pH_6.5_out) > 0) {   
  Griffith2019_pH_6.5 <- Griffith2019_pH_6.5[-Griffith2019_pH_6.5_out,] 
} 

Griffith2019_pH_6.5 <- gather(
  Griffith2019_pH_6.5,
  population,
  frequency,
  "c_a1_1" : "c_h5_1",
  factor_key = TRUE
  )

Griffith2019_pH_6.5 <- Griffith2019_pH_6.5 %>%
  select(gene, population, frequency)

Griffith2019_pH_6.5 <- Griffith2019_pH_6.5 %>%
  replace_na(value = 0)

Griffith2019_pH_6.5$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Griffith2019_pH_6.5$gene)

write_csv(Griffith2019_pH_6.5, here("data_in", "for_func", "Griffith2019_pH_6.5.csv"))

ERConstraint(
  paper = "Griffith2019",
  dataset_name = "Griffith2019_pH_6.5",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Modified LBK medium",
  generations = "1000", 
  selective_pressure = "pH 6.5",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

## Griffith2019_pH_8:
Griffith2019_pH_8 <- read_csv(
  here("data_in", "original & usable", "Griffith2019", "Griffith2019_pH_8_usable.csv")
  )

Griffith2019_pH_8 <- clean_names(Griffith2019_pH_8, case = "snake")

colnames(Griffith2019_pH_8) <- tolower(colnames(Griffith2019_pH_8))

Griffith2019_pH_8_out <- c(
  grep(out_patterns_column_gene, Griffith2019_pH_8$gene),
  grep(out_patterns_column_details, Griffith2019_pH_8$details)
  )

if (length(Griffith2019_pH_8_out) > 0) {   
  Griffith2019_pH_8 <- Griffith2019_pH_8[-Griffith2019_pH_8_out,] 
} 

Griffith2019_pH_8 <- gather(
  Griffith2019_pH_8,
  population,
  frequency,
  "c_g7_1" : "c_d11_1",
  factor_key = TRUE
  )

Griffith2019_pH_8 <- Griffith2019_pH_8 %>%
  select(gene, population, frequency)

Griffith2019_pH_8 <- Griffith2019_pH_8 %>%
  replace_na(value = 0)

Griffith2019_pH_8$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Griffith2019_pH_8$gene)

write_csv(Griffith2019_pH_8, here("data_in", "for_func", "Griffith2019_pH_8.csv"))

ERConstraint(
  paper = "Griffith2019",
  dataset_name = "Griffith2019_pH_8",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Modified LBK medium",
  generations = "1000", 
  selective_pressure = "pH 8",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
 )

# He2017
He2017 <- read_csv(
  here("data_in", "original & usable", "He2017", "He2017_usable.csv")
  )

He2017 <- clean_names(He2017, case = "snake")

colnames(He2017) <- tolower(colnames(He2017))

He2017_out <- c(
  grep(out_patterns_column_gene, He2017$gene),
  grep(out_patterns_column_details, He2017$details)
  )

if (length(He2017_out) > 0) {   
  He2017 <- He2017[-He2017_out,] 
} 

He2017 <- He2017 %>% 
  transmute(
    gene = gene, details = details, 
    b11 = rowSums(He2017[,1:2], na.rm = TRUE),
    f9 = rowSums(He2017[, 3:4], na.rm = TRUE), 
    f11 = rowSums(He2017[, 5:6], na.rm = TRUE),
    h9 = rowSums(He2017[, 7:8], na.rm = TRUE)
  )

He2017 <- gather(
  He2017, 
  population, 
  frequency, 
  "b11" : "h9",
  factor_key = TRUE
  ) 

He2017 <- He2017 %>%
  select(gene, population, frequency)

He2017 <- He2017 %>%
  replace_na(value = 0)

He2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", He2017$gene)

write_csv(He2017, here("data_in", "for_func", "He2017.csv"))

ERConstraint(
  paper = "He2017",
  dataset_name = "He2017",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Malic acid-supplemented LBK",
  generations = "2000", 
  selective_pressure = "pH 4.6 - 4.8",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Hong2011
Hong2011 <- read_csv(
  here("data_in", "original & usable", "Hong2011", "Hong2011_usable.csv")
  )

Hong2011 <- clean_names(Hong2011, case = "snake")

colnames(Hong2011) <- tolower(colnames(Hong2011))

Hong2011_out <- c(
  grep(out_patterns_column_gene, Hong2011$gene),
  grep(out_patterns_column_details, Hong2011$details)
  )

if (length(Hong2011_out) > 0) {   
  Hong2011 <- Hong2011[-Hong2011_out,] 
} 

Hong2011 <- Hong2011 %>%
  subset(details == "Coding")

Hong2011 <- Hong2011 %>%
  select(gene, population, frequency)

Hong2011 <- Hong2011 %>%
  replace_na(value = 0)

Hong2011$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Hong2011$gene)

write_csv(Hong2011, here("data_in", "for_func", "Hong2011.csv"))

ERConstraint(
  paper = "Hong2011",
  dataset_name = "Hong2011",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Galactose minimal media",
  generations = "400", 
  selective_pressure = "Galactose",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Jerison2017
Jerison2017 <- read_csv(
  here("data_in", "original & usable", "Jerison2017", "Jerison2017_usable.csv")
  )

Jerison2017 <- clean_names(Jerison2017, case = "snake")
colnames(Jerison2017) <- tolower(colnames(Jerison2017))
Jerison2017_out <- c(
  grep(out_patterns_column_gene, Jerison2017$gene),
  grep(out_patterns_column_details, Jerison2017$details)
  )

if (length(Jerison2017_out) > 0) {   
  Jerison2017 <- Jerison2017[-Jerison2017_out,] 
} 

Jerison2017 <- Jerison2017 %>%
  select(gene, population, selective_pressure, frequency)
Jerison2017 <- Jerison2017 %>%
  replace_na(value = 0)
Jerison2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Jerison2017$gene)

write_csv(Jerison2017, here("data_in", "for_func", "Jerison2017.csv"))

ERConstraint(
  paper = "Jerison2017",
  dataset_name = "Jerison2017",
  timepoint_pressure_info = "multiple",
  structure = "long",
  environment = "SC",
  generations = "500",
  selective_pressure = c("HT", "OT"), 
  species = "Sac",
  who_analyzed = "MH",
  ploidy = "haploid"
  )

# Kacar2017
Kacar2017 <- read_csv(
  here("data_in", "original & usable", "Kacar2017", "Kacar2017_usable.csv")
  )

Kacar2017 <- clean_names(Kacar2017, case = "snake")

colnames(Kacar2017) <- tolower(colnames(Kacar2017))

names(Kacar2017) <- gsub("x", "", names(Kacar2017))

Kacar2017_out <- c(
  grep(out_patterns_column_gene, Kacar2017$gene),
  grep(out_patterns_column_details, Kacar2017$details)
  )

if (length(Kacar2017_out) > 0) {   
  Kacar2017 <- Kacar2017[-Kacar2017_out,] 
} 

Kacar2017 <- Kacar2017 %>%
  select(gene, population, "500", "1000", "1500", "2000")

Kacar2017 <- Kacar2017 %>%
  replace_na(value = 0)

Kacar2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kacar2017$gene)

write_csv(Kacar2017, here("data_in", "for_func", "Kacar2017.csv"))

ERConstraint(
  paper = "Kacar2017",
  dataset_name = "Kacar2017",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "Minimal glucose medium",
  generations = c("500", "1000", "1500", "2000"), 
  selective_pressure = "Minimal glucose medium",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Native tufB replaced by inferred ancestor"
  )

# Keane2014
Keane2014 <- read_csv(
  here("data_in", "original & usable", "Keane2014", "Keane2014_usable.csv")
  )

Keane2014 = Keane2014 %>% filter(details != "intergenic")

Keane2014 <- clean_names(Keane2014, case = "snake")

colnames(Keane2014) <- tolower(colnames(Keane2014))

names(Keane2014) <- gsub("x", "", names(Keane2014))

Keane2014_out <- c(
  grep(out_patterns_column_gene, Keane2014$gene),
  grep(out_patterns_column_details, Keane2014$details)
  )

if (length(Keane2014_out) > 0) {   
  Keane2014 <- Keane2014[-Keane2014_out,] 
} 

Keane2014$`0`[Keane2014$`0` == "-"] = 0
Keane2014$`0`[Keane2014$`0` != 0] = 1
Keane2014$`440`[Keane2014$`440` != 0] = 1
Keane2014$`1100`[Keane2014$`1100` != 0] = 1
Keane2014$`1540`[Keane2014$`1540` != 0] = 1
Keane2014$`1980`[Keane2014$`1980` != 0] = 1
Keane2014$`2200`[Keane2014$`2200` != 0] = 1


Keane2014 <- Keane2014 %>%
  as.data.frame() %>%
  select(gene, population, "0","440", "1100", "1540", "1980", "2200")

Keane2014 = Keane2014 %>% filter(`0` == "0")

Keane2014$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Keane2014$gene)

write_csv(Keane2014, here("data_in", "for_func", "Keane2014.csv"))

ERConstraint(
  paper = "Keane2014",
  dataset_name = "Keane2014",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "YPD",
  generations = c("440", "1100", "1540", "1980", "2200"), 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "msh2 deleted"
  )

# Khare2015
Khare2015 <- read_csv(
  here("data_in", "original & usable", "Khare2015", "Khare2015_usable.csv")
  )

Khare2015 <- clean_names(Khare2015, case = "snake")

colnames(Khare2015) <- tolower(colnames(Khare2015))

Khare2015_out <- c(
  grep(out_patterns_column_gene, Khare2015$gene),
  grep(out_patterns_column_details, Khare2015$details)
  )

if (length(Khare2015_out) > 0) {   
  Khare2015 <- Khare2015[-Khare2015_out,] 
} 

Khare2015 <- Khare2015 %>%
  select(gene, population, frequency)

Khare2015 <- Khare2015 %>%
  replace_na(value = 0)

Khare2015$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Khare2015$gene)

write_csv(Khare2015, here("data_in", "for_func", "Khare2015.csv"))

ERConstraint(
  paper = "Khare2015",
  dataset_name = "Khare2015",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "LB",
  days = "15", 
  selective_pressure = "pvdJ spent media",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Kintses2019
Kintses2019 <- read_csv(
  here("data_in", "original & usable", "Kintses2019", "Kintses2019_usable.csv")
  )

Kintses2019 <- clean_names(Kintses2019, case = "snake")

colnames(Kintses2019) <- tolower(colnames(Kintses2019))

Kintses2019_out <- c(
  grep(out_patterns_column_gene, Kintses2019$gene),
  grep(out_patterns_column_details, Kintses2019$details)
  )

if (length(Kintses2019_out) > 0) {   
  Kintses2019 <- Kintses2019[-Kintses2019_out,] 
} 

Kintses2019 <- Kintses2019 %>%
  select(gene, population, frequency)

Kintses2019 <- Kintses2019 %>%
  replace_na(value = 0)

Kintses2019$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kintses2019$gene)

write_csv(Kintses2019, here("data_in", "for_func", "Kintses2019.csv"))

ERConstraint(
  paper = "Kintses2019",
  dataset_name = "Kintses2019",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "MS",
  generations = "120", 
  selective_pressure = "HBD-3",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Kryazhimskiy2014
## Kryazhimskiy2014_L003:
Kryazhimskiy2014_L003 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L003_usable.csv")
  )

Kryazhimskiy2014_L003 <- clean_names(Kryazhimskiy2014_L003, case = "snake")

colnames(Kryazhimskiy2014_L003) <- tolower(colnames(Kryazhimskiy2014_L003))

Kryazhimskiy2014_L003_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L003$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L003$details)
  )

if (length(Kryazhimskiy2014_L003_out) > 0) {   
  Kryazhimskiy2014_L003 <- Kryazhimskiy2014_L003[-Kryazhimskiy2014_L003_out,] 
} 
Kryazhimskiy2014_L003 <- Kryazhimskiy2014_L003 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L003 <- Kryazhimskiy2014_L003 %>%
  replace_na(value = 0)

Kryazhimskiy2014_L003$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L003$gene)

write_csv(Kryazhimskiy2014_L003, here("data_in", "for_func", "Kryazhimskiy2014_L003.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L003",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L003",
  )


## Kryazhimskiy2014_L013:
Kryazhimskiy2014_L013 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L013_usable.csv")
  )

Kryazhimskiy2014_L013 <- clean_names(Kryazhimskiy2014_L013, case = "snake")

colnames(Kryazhimskiy2014_L013) <- tolower(colnames(Kryazhimskiy2014_L013))

Kryazhimskiy2014_L013_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L013$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L013$details)
  )

if (length(Kryazhimskiy2014_L013_out) > 0) {   
  Kryazhimskiy2014_L013 <- Kryazhimskiy2014_L013[-Kryazhimskiy2014_L013_out,] 
} 

Kryazhimskiy2014_L013 <- Kryazhimskiy2014_L013 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L013 <- Kryazhimskiy2014_L013 %>%
  replace_na(value = 0)

Kryazhimskiy2014_L013$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L013$gene)

write_csv(Kryazhimskiy2014_L013, here("data_in", "for_func", "Kryazhimskiy2014_L013.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L013",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L013"
  )


## Kryazhimskiy2014_L041:
Kryazhimskiy2014_L041 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L041_usable.csv")
  )

Kryazhimskiy2014_L041 <- clean_names(Kryazhimskiy2014_L041, case = "snake")

colnames(Kryazhimskiy2014_L041) <- tolower(colnames(Kryazhimskiy2014_L041))

Kryazhimskiy2014_L041_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L041$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L041$details)
  )

if (length(Kryazhimskiy2014_L041_out) > 0) {   
  Kryazhimskiy2014_L041 <- Kryazhimskiy2014_L041[-Kryazhimskiy2014_L041_out,] 
} 

Kryazhimskiy2014_L041 <- Kryazhimskiy2014_L041 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L041 <- Kryazhimskiy2014_L041 %>%
  replace_na(value = 0)

Kryazhimskiy2014_L041$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L041$gene)

write_csv(Kryazhimskiy2014_L041, here("data_in", "for_func", "Kryazhimskiy2014_L041.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L041",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L041"
  )


## Kryazhimskiy2014_L048:
Kryazhimskiy2014_L048 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L048_usable.csv")
  )

Kryazhimskiy2014_L048 <- clean_names(Kryazhimskiy2014_L048, case = "snake")

colnames(Kryazhimskiy2014_L048) <- tolower(colnames(Kryazhimskiy2014_L048))

Kryazhimskiy2014_L048_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L048$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L048$details)
  )

if (length(Kryazhimskiy2014_L048_out) > 0) {   
  Kryazhimskiy2014_L048 <- Kryazhimskiy2014_L048[-Kryazhimskiy2014_L048_out,] 
} 

Kryazhimskiy2014_L048 <- Kryazhimskiy2014_L048 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L048 <- Kryazhimskiy2014_L048 %>%
  replace_na(value = 0)

Kryazhimskiy2014_L048$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L048$gene)

write_csv(Kryazhimskiy2014_L048, here("data_in", "for_func", "Kryazhimskiy2014_L048.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L048",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L048"
  )


## Kryazhimskiy2014_L094:
Kryazhimskiy2014_L094 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L094_usable.csv")
  )

Kryazhimskiy2014_L094 <- clean_names(Kryazhimskiy2014_L094, case = "snake")

colnames(Kryazhimskiy2014_L094) <- tolower(colnames(Kryazhimskiy2014_L094))

Kryazhimskiy2014_L094_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L094$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L094$details)
  )

if (length(Kryazhimskiy2014_L094_out) > 0) {   
  Kryazhimskiy2014_L094 <- Kryazhimskiy2014_L094[-Kryazhimskiy2014_L094_out,] 
} 

Kryazhimskiy2014_L094 <- Kryazhimskiy2014_L094 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L094 <- Kryazhimskiy2014_L094 %>%
  replace_na(value = 0)

Kryazhimskiy2014_L094$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L094$gene)

write_csv(Kryazhimskiy2014_L094, here("data_in", "for_func", "Kryazhimskiy2014_L094.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L094",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L094"
  )


## Kryazhimskiy2014_L096a:
Kryazhimskiy2014_L096a <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L096a_usable.csv")
  )

Kryazhimskiy2014_L096a <- clean_names(Kryazhimskiy2014_L096a, case = "snake")

colnames(Kryazhimskiy2014_L096a) <- tolower(colnames(Kryazhimskiy2014_L096a))

Kryazhimskiy2014_L096a_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L096a$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L096a$details)
  )

if (length(Kryazhimskiy2014_L096a_out) > 0) {   
  Kryazhimskiy2014_L096a <- Kryazhimskiy2014_L096a[-Kryazhimskiy2014_L096a_out,] 
} 

Kryazhimskiy2014_L096a <- Kryazhimskiy2014_L096a %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L096a <- Kryazhimskiy2014_L096a %>%
  replace_na(value = 0)

Kryazhimskiy2014_L096a$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L096a$gene)

write_csv(Kryazhimskiy2014_L096a, here("data_in", "for_func", "Kryazhimskiy2014_L096a.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L096a",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L096a"
  )


## Kryazhimskiy2014_L096b:
Kryazhimskiy2014_L096b <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L096b_usable.csv")
  )

Kryazhimskiy2014_L096b <- clean_names(Kryazhimskiy2014_L096b, case = "snake")

colnames(Kryazhimskiy2014_L096b) <- tolower(colnames(Kryazhimskiy2014_L096b))

Kryazhimskiy2014_L096b_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L096b$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L096b$details)
  )

if (length(Kryazhimskiy2014_L096b_out) > 0) {   
  Kryazhimskiy2014_L096b <- Kryazhimskiy2014_L096b[-Kryazhimskiy2014_L096b_out,] 
} 

Kryazhimskiy2014_L096b <- Kryazhimskiy2014_L096b %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L096b <- Kryazhimskiy2014_L096b %>%
  replace_na(value = 0)

Kryazhimskiy2014_L096b$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L096b$gene)

write_csv(Kryazhimskiy2014_L096b, here("data_in", "for_func", "Kryazhimskiy2014_L096b.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L096b",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L096b"
  )


## Kryazhimskiy2014_L098:
Kryazhimskiy2014_L098 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L098_usable.csv")
  )

Kryazhimskiy2014_L098 <- clean_names(Kryazhimskiy2014_L098, case = "snake")

colnames(Kryazhimskiy2014_L098) <- tolower(colnames(Kryazhimskiy2014_L098))

Kryazhimskiy2014_L098_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L098$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L098$details)
  )

if (length(Kryazhimskiy2014_L098_out) > 0) {   
  Kryazhimskiy2014_L098 <- Kryazhimskiy2014_L098[-Kryazhimskiy2014_L098_out,] 
} 

Kryazhimskiy2014_L098 <- Kryazhimskiy2014_L098 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L098 <- Kryazhimskiy2014_L098 %>%
  replace_na(value = 0)

Kryazhimskiy2014_L098$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L098$gene)

write_csv(Kryazhimskiy2014_L098, here("data_in", "for_func", "Kryazhimskiy2014_L098.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L098",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
 selective_pressure = "YPD",
 species = "Sac",
 who_analyzed = "TL",
 ploidy = "haploid",
 strain_info = "Founder L098"
 )


## Kryazhimskiy2014_L102:
Kryazhimskiy2014_L102 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L102_usable.csv")
  )

Kryazhimskiy2014_L102 <- clean_names(Kryazhimskiy2014_L102, case = "snake")

colnames(Kryazhimskiy2014_L102) <- tolower(colnames(Kryazhimskiy2014_L102))

Kryazhimskiy2014_L102_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L102$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L102$details)
  )

if (length(Kryazhimskiy2014_L102_out) > 0) {   
  Kryazhimskiy2014_L102 <- Kryazhimskiy2014_L102[-Kryazhimskiy2014_L102_out,] 
} 

Kryazhimskiy2014_L102 <- Kryazhimskiy2014_L102 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L102 <- Kryazhimskiy2014_L102 %>%
  replace_na(value = 0)

Kryazhimskiy2014_L102$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L102$gene)

write_csv(Kryazhimskiy2014_L102, here("data_in", "for_func", "Kryazhimskiy2014_L102.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L102",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L102"
  )


## Kryazhimskiy2014_L102a:
Kryazhimskiy2014_L102a <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_L102a_usable.csv")
  )

Kryazhimskiy2014_L102a <- clean_names(Kryazhimskiy2014_L102a, case = "snake")

colnames(Kryazhimskiy2014_L102a) <- tolower(colnames(Kryazhimskiy2014_L102a))

Kryazhimskiy2014_L102a_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_L102a$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_L102a$details)
  )

if (length(Kryazhimskiy2014_L102a_out) > 0) {   
  Kryazhimskiy2014_L102a <- Kryazhimskiy2014_L102a[-Kryazhimskiy2014_L102a_out,] 
} 

Kryazhimskiy2014_L102a <- Kryazhimskiy2014_L102a %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_L102a <- Kryazhimskiy2014_L102a %>%
  replace_na(value = 0)

Kryazhimskiy2014_L102a$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_L102a$gene)

write_csv(Kryazhimskiy2014_L102a, here("data_in", "for_func", "Kryazhimskiy2014_L102a.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_L102a",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder L102a"
  )


## Kryazhimskiy2014_S002:
Kryazhimskiy2014_S002 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_S002_usable.csv")
  )

Kryazhimskiy2014_S002 <- clean_names(Kryazhimskiy2014_S002, case = "snake")

colnames(Kryazhimskiy2014_S002) <- tolower(colnames(Kryazhimskiy2014_S002))

Kryazhimskiy2014_S002_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_S002$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_S002$details)
  )

if (length(Kryazhimskiy2014_S002_out) > 0) {   
  Kryazhimskiy2014_S002 <- Kryazhimskiy2014_S002[-Kryazhimskiy2014_S002_out,] 
} 

Kryazhimskiy2014_S002 <- Kryazhimskiy2014_S002 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_S002 <- Kryazhimskiy2014_S002 %>%
  replace_na(value = 0)

Kryazhimskiy2014_S002$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_S002$gene)

write_csv(Kryazhimskiy2014_S002, here("data_in", "for_func", "Kryazhimskiy2014_S002.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_S002",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder S002"
  )


## Kryazhimskiy2014_S028:
Kryazhimskiy2014_S028 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_S028_usable.csv")
  )

Kryazhimskiy2014_S028 <- clean_names(Kryazhimskiy2014_S028, case = "snake")

colnames(Kryazhimskiy2014_S028) <- tolower(colnames(Kryazhimskiy2014_S028))

Kryazhimskiy2014_S028_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_S028$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_S028$details)
  )

if (length(Kryazhimskiy2014_S028_out) > 0) {   
  Kryazhimskiy2014_S028 <- Kryazhimskiy2014_S028[-Kryazhimskiy2014_S028_out,] 
} 

Kryazhimskiy2014_S028 <- Kryazhimskiy2014_S028 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_S028 <- Kryazhimskiy2014_S028 %>%
  replace_na(value = 0)

Kryazhimskiy2014_S028$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_S028$gene)

write_csv(Kryazhimskiy2014_S028, here("data_in", "for_func", "Kryazhimskiy2014_S028.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_S028",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder S028"
  )


## Kryazhimskiy2014_S121:
Kryazhimskiy2014_S121 <- read_csv(
  here("data_in", "original & usable", "Kryazhimskiy2014", "Kryazhimskiy2014_S121_usable.csv")
  )

Kryazhimskiy2014_S121 <- clean_names(Kryazhimskiy2014_S121, case = "snake")

colnames(Kryazhimskiy2014_S121) <- tolower(colnames(Kryazhimskiy2014_S121))

Kryazhimskiy2014_S121_out <- c(
  grep(out_patterns_column_gene, Kryazhimskiy2014_S121$gene),
  grep(out_patterns_column_details, Kryazhimskiy2014_S121$details)
  )

if (length(Kryazhimskiy2014_S121_out) > 0) {   
  Kryazhimskiy2014_S121 <- Kryazhimskiy2014_S121[-Kryazhimskiy2014_S121_out,] 
} 

Kryazhimskiy2014_S121 <- Kryazhimskiy2014_S121 %>%
  subset(is.na(distance_to_gene)) %>%
  select(gene, population, frequency)

Kryazhimskiy2014_S121 <- Kryazhimskiy2014_S121 %>%
  replace_na(value = 0)

Kryazhimskiy2014_S121$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Kryazhimskiy2014_S121$gene)

write_csv(Kryazhimskiy2014_S121, here("data_in", "for_func", "Kryazhimskiy2014_S121.csv"))

ERConstraint(
  paper = "Kryazhimskiy2014",
  dataset_name = "Kryazhimskiy2014_S121",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  generations = "500", 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Founder S121"
  )

# KuzdzalFick2018
KuzdzalFick2018 <- read_csv(
  here("data_in", "original & usable", "KuzdzalFick2018", "KuzdzalFick2018_usable.csv")
  )

KuzdzalFick2018 <- clean_names(KuzdzalFick2018, case = "snake")

colnames(KuzdzalFick2018) <- tolower(colnames(KuzdzalFick2018))

KuzdzalFick2018_out <- c(
  grep(out_patterns_column_gene, KuzdzalFick2018$gene),
  grep(out_patterns_column_details, KuzdzalFick2018$details)
  )

if (length(KuzdzalFick2018_out) > 0) {   
  KuzdzalFick2018 <- KuzdzalFick2018[-KuzdzalFick2018_out,] 
} 

KuzdzalFick2018 <- KuzdzalFick2018 %>%
  select(gene, population, frequency)

KuzdzalFick2018 <- KuzdzalFick2018 %>%
  replace_na(value = 0)

KuzdzalFick2018$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", KuzdzalFick2018$gene)

write_csv(KuzdzalFick2018, here("data_in", "for_func", "KuzdzalFick2018.csv"))

ERConstraint(
  paper = "KuzdzalFick2018",
  dataset_name = "KuzdzalFick2018",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YPD",
  days = "28", 
  selective_pressure = "M9 minimal medium",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "TBR1, a segregant obtained by multiple crosses of baking strains."
  )


# Lang2013
Lang2013 <- read_csv(
  here("data_in", "original & usable", "Lang2013", "Lang2013_usable.csv")
  )

Lang2013 <- clean_names(Lang2013, case = "snake")

colnames(Lang2013) <- tolower(colnames(Lang2013))

names(Lang2013) <- gsub("p1_", "", names(Lang2013))

Lang2013_out <- c(
  grep(",", Lang2013$gene),
  grep("/", Lang2013$gene),
  grep("\\[", Lang2013$gene),
  grep("prophage", Lang2013$gene), 
  grep(out_patterns_column_details,
  Lang2013$details)
  )

if (length(Lang2013_out) > 0) {
  Lang2013 <- Lang2013[-Lang2013_out,]
}

Lang2013 <- Lang2013 %>%
  select(gene,
         population,
         '140',
         '240',
         '335',
         '415',
         '505',
         '585',
         '665',
         '745',
         '825',
         '910',
         '1000'
         )

Lang2013 <- Lang2013 %>%
  replace_na(value = 0)

Lang2013$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Lang2013$gene)

write_csv(Lang2013, here("data_in", "for_func", "Lang2013.csv"))

ERConstraint(
  paper = "Lang2013",
  dataset_name = "Lang2013",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "YPD",
  generations = c('140', '240', '335', '415', '505', '585', '665', '745', '825', '910', '1000'), 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "MH",
  ploidy = "haploid"
  )

# Long2017
Long2017 <- read_csv(
  here("data_in", "original & usable", "Long2017", "Long2017_usable.csv")
  )

Long2017 <- clean_names(Long2017, case = "snake")

colnames(Long2017) <- tolower(colnames(Long2017))

Long2017_out <- c(
  grep(out_patterns_column_gene, Long2017$gene),
  grep(out_patterns_column_details, Long2017$details)
  )

if (length(Long2017_out) > 0) {
  Long2017 <- Long2017[-Long2017_out,] 
} 

Long2017 <- gather(Long2017, population, frequency, "ale1":"ale10", factor_key=TRUE)

Long2017 <- Long2017 %>%
  select(gene, population, frequency)

Long2017 <- Long2017 %>%
  replace_na(value = 0)

Long2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Long2017$gene)

write_csv(Long2017, here("data_in", "for_func", "Long2017.csv"))

ERConstraint(
  paper = "Long2017",
  dataset_name = "Long2017",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Glucose minimal media",
  days = "50",
  selective_pressure = "Glucose minimal media", 
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid",
  strain_info = "pgi deleted"
  )

# McCloskey2018
McCloskey2018 <- read_csv(here("data_in", "original & usable", "McCloskey2018", "McCloskey2018_usable.csv"))
McCloskey2018 <- clean_names(McCloskey2018, case = "snake")
colnames(McCloskey2018) <- tolower(colnames(McCloskey2018))

McCloskey2018_out <- c(
  grep(out_patterns_column_gene, McCloskey2018$gene),
  grep(out_patterns_column_details, McCloskey2018$details)
  )

if (length(McCloskey2018_out) > 0) {   
  McCloskey2018 <- McCloskey2018[-McCloskey2018_out,] 
}

McCloskey2018 <- McCloskey2018 %>%
  select(gene, population, frequency)

McCloskey2018 <- McCloskey2018 %>%
  replace_na(value = 0)

McCloskey2018$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", McCloskey2018$gene)

## McCloskey2018_gnd:
McCloskey2018_gnd_rows <- grep("gnd", McCloskey2018$population)

McCloskey2018_gnd <- McCloskey2018[McCloskey2018_gnd_rows,]

write_csv(McCloskey2018_gnd, here("data_in", "for_func", "McCloskey2018_gnd.csv"))

ERConstraint(
  paper = "McCloskey2018",
  dataset_name = "McCloskey2018_gnd",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  days = "35", 
  selective_pressure = "M9 minimal medium",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "gnd"
  )


## McCloskey2018_pgi:
McCloskey2018_pgi_rows <- grep("pgi", McCloskey2018$population)
McCloskey2018_pgi <- McCloskey2018[McCloskey2018_pgi_rows,]

write_csv(McCloskey2018_pgi, here("data_in", "for_func", "McCloskey2018_pgi.csv"))

ERConstraint(
  paper = "McCloskey2018",
  dataset_name = "McCloskey2018_pgi",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  days = "35", 
  selective_pressure = "M9 minimal medium",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "pgi"
  )

## McCloskey2018_ptsHIcrr:
McCloskey2018_ptsHIcrr_rows <- grep("ptsHIcrr", McCloskey2018$population)

McCloskey2018_ptsHIcrr <- McCloskey2018[McCloskey2018_ptsHIcrr_rows,]

write_csv(McCloskey2018_ptsHIcrr, here("data_in", "for_func", "McCloskey2018_ptsHIcrr.csv"))

ERConstraint(
  paper = "McCloskey2018",
  dataset_name = "McCloskey2018_ptsHIcrr",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  days = "35", 
  selective_pressure = "M9 minimal medium",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info ="ptsHIcrr"
  )

##  McCloskey2018_sdhCB:
McCloskey2018_sdhCB_rows <- grep("sdhCB", McCloskey2018$population)

McCloskey2018_sdhCB <- McCloskey2018[McCloskey2018_sdhCB_rows,]

write_csv(McCloskey2018_sdhCB, here("data_in", "for_func", "McCloskey2018_sdhCB.csv"))

ERConstraint(
  paper = "McCloskey2018",
  dataset_name = "McCloskey2018_sdhCB",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  days = "35", 
  selective_pressure = "M9 minimal medium",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "sdhCB"
  )


## McCloskey2018_tpiAE:
McCloskey2018_tpiAE_rows <- grep("tpiAE", McCloskey2018$population)

McCloskey2018_tpiAE <- McCloskey2018[McCloskey2018_tpiAE_rows,]

write_csv(McCloskey2018_tpiAE, here("data_in", "for_func", "McCloskey2018_tpiAE.csv"))

ERConstraint(
  paper = "McCloskey2018",
  dataset_name = "McCloskey2018_tpiAE",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  days = "35", 
  selective_pressure = "M9 minimal medium",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "tpiAE"
  )

## McCloskey2018_Evo:
McCloskey2018_Evo_rows <- grep("Evo", McCloskey2018$population)

McCloskey2018_Evo <- McCloskey2018[McCloskey2018_Evo_rows,]

write_csv(McCloskey2018_Evo, here("data_in", "for_func", "McCloskey2018_Evo.csv"))

ERConstraint(
  paper = "McCloskey2018",
  dataset_name = "McCloskey2018_Evo",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  days = "35", 
  selective_pressure = "M9 minimal medium",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid",
  strain_info = "Evo"
  )
# Monk2016
Monk2016 <- read_csv(
  here("data_in", "original & usable", "Monk2016", "Monk2016_usable.csv")
  )

Monk2016 <- clean_names(Monk2016, case = "snake")

colnames(Monk2016) <- tolower(colnames(Monk2016))

Monk2016_out <- c(
  grep(out_patterns_column_gene, Monk2016$gene),
  grep(out_patterns_column_details, Monk2016$details)
  )

if (length(Monk2016_out) > 0) {   
  Monk2016 <- Monk2016[-Monk2016_out,] 
} 

Monk2016 <- Monk2016 %>%
  select(gene, population, frequency)

Monk2016 <- Monk2016 %>%
  replace_na(value = 0)

Monk2016$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Monk2016$gene)

write_csv(Monk2016, here("data_in", "for_func", "Monk2016.csv"))

ERConstraint(
  paper = "Monk2016",
  dataset_name = "Monk2016",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "LB alternating with M9GBT (M9 + glucose + biotin + thiamine)",
  days = "18", 
  selective_pressure = "TAG amber codons replaced with TAA stop codons",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Morgenthaler2019
Morgenthaler2019 <- read_csv(
  here("data_in", "original & usable", "Morgenthaler2019", "Morgenthaler2019_usable.csv")
  )

Morgenthaler2019 <- clean_names(Morgenthaler2019, case = "snake")

colnames(Morgenthaler2019) <- tolower(colnames(Morgenthaler2019))

Morgenthaler2019_out <- c(
  grep(out_patterns_column_gene, Morgenthaler2019$gene),
  grep(out_patterns_column_details, Morgenthaler2019$details)
  )

if (length(Morgenthaler2019_out) > 0) {   
  Morgenthaler2019 <- Morgenthaler2019[-Morgenthaler2019_out,] 
}

Morgenthaler2019 <- Morgenthaler2019 %>%
  transmute(gene, population, "42" = "day_42", "50" = "day_50")

Morgenthaler2019 <- Morgenthaler2019 %>%
  replace_na(value = 0)

Morgenthaler2019$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Morgenthaler2019$gene)


write_csv(Morgenthaler2019, here("data_in", "for_func", "Morgenthaler2019.csv"))

ERConstraint(
  paper = "Morgenthaler2019",
  dataset_name = "Morgenthaler2019",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "M9",
  days = c("42", "50"),
  selective_pressure = "M9", 
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Mundhada2017
Mundhada2017 <- read_csv(
  here("data_in", "original & usable", "Mundhada2017", "Mundhada2017_usable.csv")
  )

Mundhada2017 <- clean_names(Mundhada2017, case = "snake")

colnames(Mundhada2017) <- tolower(colnames(Mundhada2017))

Mundhada2017_out <- c(
  grep(out_patterns_column_gene, Mundhada2017$gene),
  grep(out_patterns_column_details, Mundhada2017$details)
  )

if (length(Mundhada2017_out) > 0) {   
  Mundhada2017 <- Mundhada2017[-Mundhada2017_out,] 
} 

Mundhada2017 <- Mundhada2017 %>% 
  transmute(
    gene = gene, 
    details = details, 
    a3 = rowSums(Mundhada2017[, 3:4], na.rm = TRUE),
    a4 = rowSums(Mundhada2017[, 5:6], na.rm = TRUE), 
    a5 = rowSums(Mundhada2017[, 7:8], na.rm = TRUE)
  ) %>%
  gather(population, frequency, "a3" : "a5", factor_key=TRUE) %>%
  select(gene, population, frequency)

Mundhada2017 <- Mundhada2017 %>%
  replace_na(value = 0)

Mundhada2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Mundhada2017$gene)

write_csv(Mundhada2017, here("data_in", "for_func", "Mundhada2017.csv"))

ERConstraint(
  paper = "Mundhada2017",
  dataset_name = "Mundhada2017",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium + glucose + acetate",
  generations = "650", 
  selective_pressure = "glucose + acetate",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Payen2016
Payen2016 <- read_csv(here("data_in", "original & usable", "Payen2016", "Payen2016_usable.csv"))

Payen2016 <- clean_names(Payen2016, case = "snake")

colnames(Payen2016) <- tolower(colnames(Payen2016))

Payen2016_out <- c(
  grep(out_patterns_column_gene, Payen2016$gene),
  grep(out_patterns_column_details, Payen2016$details)
  )

if (length(Payen2016_out) > 0) {   
  Payen2016 <- Payen2016[-Payen2016_out,] 
}

Payen2016 <- Payen2016 %>%  
  filter(frequency != "clone", frequency != "0")

Payen2016 <- Payen2016 %>%
  replace_na(value = 0)

Payen2016$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Payen2016$gene)

## Payen2016_dip:
Payen2016_dip <- Payen2016 %>%
  filter(ploidy == "diploid") %>%
  select(gene, population, selective_pressure, frequency)

write_csv(Payen2016_dip, here("data_in", "for_func", "Payen2016_dip.csv"))

ERConstraint(
  paper = "Payen2016",
  dataset_name = "Payen2016_dip",
  timepoint_pressure_info = "multiple",
  structure = "long",
  environment = "phosphate",
  generations = "20",
  selective_pressure = c("phosphate"),
  species = "Sac", 
  who_analyzed = "MH",
  ploidy = "diploid",
  strain_info = "diploid"
  )


## Payen2016_hap:
Payen2016_hap <- Payen2016 %>%
  filter(ploidy == "haploid") %>%
  select(gene, population, selective_pressure, frequency)

write_csv(Payen2016_hap, here("data_in", "for_func", "Payen2016_hap.csv"))

ERConstraint(
  paper = "Payen2016",
  dataset_name = "Payen2016_hap",
  timepoint_pressure_info = "multiple",
  structure = "long",
  environment = "LB",
  generations = "20",
  selective_pressure = c("phosphate"),
  species = "Sac", 
  who_analyzed = "MH",
  ploidy = "haploid",
  strain_info = "haploid"
  )

# Sandberg2014
Sandberg2014 <- read_csv(
  here("data_in", "original & usable", "Sandberg2014", "Sandberg2014_usable.csv")
  )
Sandberg2014 <- clean_names(Sandberg2014, case = "snake")

colnames(Sandberg2014) <- tolower(colnames(Sandberg2014))

Sandberg2014_out <- c(
  grep(out_patterns_column_gene, Sandberg2014$gene),
  grep(out_patterns_column_details, Sandberg2014$details)
  )

if (length(Sandberg2014_out) > 0) {
  Sandberg2014 <- Sandberg2014[-Sandberg2014_out,] 
} 

Sandberg2014 <- gather(
  Sandberg2014,
  population,
  frequency,
  "ale_number_1":"ale_number_10",
  factor_key = TRUE
  )

Sandberg2014 <- Sandberg2014 %>%
  select(gene, population, frequency)

Sandberg2014 <- Sandberg2014 %>%
  replace_na(value = 0)

Sandberg2014$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandberg2014$gene)

write_csv(Sandberg2014, here("data_in", "for_func", "Sandberg2014.csv"))

ERConstraint(
  paper = "Sandberg2014",
  dataset_name = "Sandberg2014",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Glucose minimal media",
  days = "45", 
  selective_pressure = "heat (42 degrees C)",
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid"
  )

# Sandberg2016
Sandberg2016 <- read_csv(
  here("data_in", "original & usable", "Sandberg2016", "Sandberg2016_usable.csv")
  )

Sandberg2016 <- clean_names(Sandberg2016, case = "snake")

colnames(Sandberg2016) <- tolower(colnames(Sandberg2016))

Sandberg2016_out <- c(
  grep(out_patterns_column_gene, Sandberg2016$gene),
  grep(out_patterns_column_details, Sandberg2016$details)
  )

if (length(Sandberg2016_out) > 0) {   
  Sandberg2016 <- Sandberg2016[-Sandberg2016_out,] 
} 

Sandberg2016 <- Sandberg2016 %>%
  transmute(
    gene,
    population,
    "23" = Sandberg2016$"flask_23",
    "58" = Sandberg2016$"flask_58",
    "133" = Sandberg2016$"flask_133"
    )

Sandberg2016 <- Sandberg2016 %>%
  replace_na(value = 0)

Sandberg2016$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandberg2016$gene)

write_csv(Sandberg2016, here("data_in", "for_func", "Sandberg2016.csv"))

ERConstraint(
  paper = "Sandberg2016",
  dataset_name = "Sandberg2016",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "Davis minimal medium",
  selective_pressure = "13C isotope", 
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid",
  flasks = c("23", "58", "133")
  )

# Sandberg2017  
## Sandberg2017_ac:
Sandberg2017_ac <- read_csv(
  here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_ac_usable.csv")
  )

Sandberg2017_ac <- clean_names(Sandberg2017_ac, case = "snake")

colnames(Sandberg2017_ac) <- tolower(colnames(Sandberg2017_ac))

Sandberg2017_ac_out <- c(
  grep(out_patterns_column_gene, Sandberg2017_ac$gene),
  grep(out_patterns_column_details, Sandberg2017_ac$details)
  )

if (length(Sandberg2017_ac_out) > 0) {   
  Sandberg2017_ac <- Sandberg2017_ac[-Sandberg2017_ac_out,] 
} 

Sandberg2017_ac <- gather(
  Sandberg2017_ac, 
  population, 
  frequency, 
  "a1_f56_i1_r1" : "a4_f55_i1_r1",
  factor_key = TRUE
  )

Sandberg2017_ac <- Sandberg2017_ac %>%
  select(gene, population, frequency)

Sandberg2017_ac <- Sandberg2017_ac %>%
  replace_na(value = 0)

Sandberg2017_ac$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandberg2017_ac$gene)

write_csv(Sandberg2017_ac, here("data_in", "for_func", "Sandberg2017_ac.csv"))

ERConstraint(
  paper = "Sandberg2017",
  dataset_name = "Sandberg2017_ac",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium + acetate",
  generations = "1000", 
  selective_pressure = "acetate",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

## Sandberg2017_glu_ac:
Sandberg2017_glu_ac <- read_csv(
  here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_ac_usable.csv")
  )

Sandberg2017_glu_ac <- clean_names(Sandberg2017_glu_ac, case = "snake")

colnames(Sandberg2017_glu_ac) <- tolower(colnames(Sandberg2017_glu_ac))

Sandberg2017_glu_ac_out <- c(
  grep(out_patterns_column_gene, Sandberg2017_glu_ac$gene),
  grep(out_patterns_column_details, Sandberg2017_glu_ac$details)
  )

if (length(Sandberg2017_glu_ac_out) > 0) {   
  Sandberg2017_glu_ac <- Sandberg2017_glu_ac[-Sandberg2017_glu_ac_out,] 
} 

Sandberg2017_glu_ac <- Sandberg2017_glu_ac %>% 
  transmute(
    gene = gene, details = details, 
    a7 = rowSums(Sandberg2017_glu_ac[, 4:15], na.rm = TRUE),
    a8 = rowSums(Sandberg2017_glu_ac[, 16:28], na.rm = TRUE), 
    a9 = rowSums(Sandberg2017_glu_ac[, 29:38], na.rm = TRUE)
  )

Sandberg2017_glu_ac <- gather(
  Sandberg2017_glu_ac,
  population, frequency,
  "a7" : "a9",
  factor_key = TRUE 
  )

Sandberg2017_glu_ac <- Sandberg2017_glu_ac %>%
  select(gene, population, frequency)

Sandberg2017_glu_ac <- Sandberg2017_glu_ac %>%
  replace_na(value = 0)

Sandberg2017_glu_ac$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandberg2017_glu_ac$gene)

write_csv(Sandberg2017_glu_ac, here("data_in", "for_func", "Sandberg2017_glu_ac.csv"))

ERConstraint(
  paper = "Sandberg2017",
  dataset_name = "Sandberg2017_glu_ac",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium + glucose + acetate",
  generations = "650", 
  selective_pressure = "glucose + acetate",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

## Sandberg2017_glu_gly:
Sandberg2017_glu_gly <- read_csv(
  here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_gly_usable.csv")
  )

Sandberg2017_glu_gly <- clean_names(Sandberg2017_glu_gly, case = "snake")

colnames(Sandberg2017_glu_gly) <- tolower(colnames(Sandberg2017_glu_gly))

Sandberg2017_glu_gly_out <- c(
  grep(out_patterns_column_gene, Sandberg2017_glu_gly$gene),
  grep(out_patterns_column_details, Sandberg2017_glu_gly$details)
  )

if (length(Sandberg2017_glu_gly_out) > 0) {   
  Sandberg2017_glu_gly <- Sandberg2017_glu_gly[-Sandberg2017_glu_gly_out,] 
} 

Sandberg2017_glu_gly <- Sandberg2017_glu_gly %>% 
  transmute(
    gene = gene, details = details, 
    a4 = rowSums(Sandberg2017_glu_gly[, 4:7], na.rm = TRUE),
    a5 = rowSums(Sandberg2017_glu_gly[, 8:12], na.rm = TRUE), 
    a6 = rowSums(Sandberg2017_glu_gly[, 13:16], na.rm = TRUE)
)

Sandberg2017_glu_gly <- gather(
  Sandberg2017_glu_gly,
  population,
  frequency,
  "a4" : "a6",
  factor_key = TRUE
  )

Sandberg2017_glu_gly <- Sandberg2017_glu_gly %>%
  select(gene, population, frequency)

Sandberg2017_glu_gly <- Sandberg2017_glu_gly %>%
  replace_na(value = 0)

Sandberg2017_glu_gly$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandberg2017_glu_gly$gene)

write_csv(Sandberg2017_glu_gly, here("data_in", "for_func", "Sandberg2017_glu_gly.csv"))

ERConstraint(
  paper = "Sandberg2017",
  dataset_name = "Sandberg2017_glu_gly",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium + glucose + glycerol",
  generations = "1170", 
  selective_pressure = "glucose + glycerol",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

## Sandberg2017_glu_xyl:
Sandberg2017_glu_xyl <- read_csv(
  here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_glu_xyl_usable.csv")
  )

Sandberg2017_glu_xyl <- clean_names(Sandberg2017_glu_xyl, case = "snake")

colnames(Sandberg2017_glu_xyl) <- tolower(colnames(Sandberg2017_glu_xyl))

Sandberg2017_glu_xyl_out <- c(
  grep(out_patterns_column_gene, Sandberg2017_glu_xyl$gene),
  grep(out_patterns_column_details, Sandberg2017_glu_xyl$details)
  )

if (length(Sandberg2017_glu_xyl_out) > 0) {   
  Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl[-Sandberg2017_glu_xyl_out,] 
} 

Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl %>% 
  transmute(
    gene = gene, details = details, 
    a0 = rowSums(Sandberg2017_glu_xyl[, 3], na.rm = TRUE),
    a1 = rowSums(Sandberg2017_glu_xyl[, 4:13], na.rm = TRUE), 
    a2 = rowSums(Sandberg2017_glu_xyl[, 14:19], na.rm = TRUE),
    a3 = rowSums(Sandberg2017_glu_xyl[, 20:ncol(Sandberg2017_glu_xyl)], na.rm = TRUE)
  )

Sandberg2017_glu_xyl <- gather(
  Sandberg2017_glu_xyl,
  population,
  frequency,
  "a1" : "a3",
  factor_key = TRUE
  )

Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl %>%
  select(gene, population, frequency)

Sandberg2017_glu_xyl <- Sandberg2017_glu_xyl %>%
  replace_na(value = 0)

Sandberg2017_glu_xyl$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandberg2017_glu_xyl$gene)

write_csv(Sandberg2017_glu_xyl, here("data_in", "for_func", "Sandberg2017_glu_xyl.csv"))

ERConstraint(
  paper = "Sandberg2017",
  dataset_name = "Sandberg2017_glu_xyl",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium + glucose + xylose",
  generations = "1180", 
  selective_pressure = "glucose + xylose",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

## Sandberg2017_xyl:
Sandberg2017_xyl <- read_csv(
  here("data_in", "original & usable", "Sandberg2017", "Sandberg2017_xyl_usable.csv")
  )

Sandberg2017_xyl <- clean_names(Sandberg2017_xyl, case = "snake")

colnames(Sandberg2017_xyl) <- tolower(colnames(Sandberg2017_xyl))

Sandberg2017_xyl_out <- c(
  grep(out_patterns_column_gene, Sandberg2017_xyl$gene),
  grep(out_patterns_column_details, Sandberg2017_xyl$details)
  )

if (length(Sandberg2017_xyl_out) > 0) {   
  Sandberg2017_xyl <- Sandberg2017_xyl[-Sandberg2017_xyl_out,] 
} 

Sandberg2017_xyl <- gather(
  Sandberg2017_xyl,
  population,
  frequency,
  "a1_f116_i1_r1" : "a4_f113_i1_r1",
  factor_key = TRUE
  )

Sandberg2017_xyl <- Sandberg2017_xyl %>%
  select(gene, population, frequency)

Sandberg2017_xyl <- Sandberg2017_xyl %>%
  replace_na(value = 0)

Sandberg2017_xyl$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandberg2017_xyl$gene)

write_csv(Sandberg2017_xyl, here("data_in", "for_func", "Sandberg2017_xyl.csv"))

ERConstraint(
  paper = "Sandberg2017",
  dataset_name = "Sandberg2017_xyl",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium + xylose",
  generations = "1000", 
  selective_pressure = "xylose",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Sandra2016
Sandra2016 <- read_csv(
  here("data_in", "original & usable", "Sandra2016", "Sandra2016_usable.csv")
  )

Sandra2016 <- clean_names(Sandra2016, case = "snake")

colnames(Sandra2016) <- tolower(colnames(Sandra2016))

Sandra2016_out <- c(
  grep(out_patterns_column_gene, Sandra2016$gene),
  grep(out_patterns_column_details, Sandra2016$details)
  )

if (length(Sandra2016_out) > 0) {   
  Sandra2016 <- Sandra2016[-Sandra2016_out,] 
} 

Sandra2016 <- Sandra2016 %>%
  select(gene, population, frequency)
Sandra2016 <- Sandra2016 %>%
  replace_na(value = 0)
Sandra2016$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sandra2016$gene)

write_csv(Sandra2016, here("data_in", "for_func", "Sandra2016.csv"))
ERConstraint(
  paper = "Sandra2016",
  dataset_name = "Sandra2016",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  days = "8", 
  selective_pressure = "Combined UVA and UVB irradiation",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Saxer
## Saxer2014_ancestorA:
Saxer2014_ancestorA <- read_csv(
  here("data_in", "original & usable", "Saxer2014", "Saxer2014_ancestorA_usable.csv")
  )

Saxer2014_ancestorA <- clean_names(Saxer2014, case = "snake")

colnames(Saxer2014_ancestorA) <- tolower(colnames(Saxer2014_ancestorA))

Saxer2014_ancestorA_out <- c(
  grep(out_patterns_column_gene, Saxer2014_ancestorA$gene), 
  grep(out_patterns_column_details, Saxer2014_ancestorA$details)
  )

if (length(Saxer2014_ancestorA_out) > 0) {   
  Saxer2014_ancestorA <- Saxer2014[-Saxer2014_ancestorA_out,] 
} 

Saxer2014_ancestorA <- Saxer2014_ancestorA %>%
  select(gene, population, frequency)

Saxer2014_ancestorA <- Saxer2014_ancestorA %>%
  replace_na(value = 0)

Saxer2014_ancestorA$frequency <- as.integer(as.logical(Saxer2014_ancestorA$frequency))

Saxer2014_ancestorA$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Saxer2014_ancestorA$gene)

write_csv(Saxer2014_ancestorA, here("data_in", "for_func", "Saxer2014_ancestorA.csv"))

ERConstraint(
  paper = "Saxer2014",
  dataset_name = "Saxer2014_ancestorA",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "BHI",
  generations = "765", 
  selective_pressure = "Rich media (BBL BHI)",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

## Saxer2014_ancestorB:
Saxer2014_ancestorB <- read_csv(
  here("data_in", "original & usable", "Saxer2014", "Saxer2014_ancestorB_usable.csv")
  )

Saxer2014_ancestorB <- clean_names(Saxer2014, case = "snake")

colnames(Saxer2014_ancestorB) <- tolower(colnames(Saxer2014_ancestorB))

Saxer2014_ancestorB_out <- c(
  grep(out_patterns_column_gene, Saxer2014_ancestorB$gene),
  grep(out_patterns_column_details, Saxer2014_ancestorB$details)
  )

if (length(Saxer2014_ancestorB_out) > 0) {   
  Saxer2014_ancestorB <- Saxer2014[-Saxer2014_ancestorB_out,] 
} 

Saxer2014_ancestorB <- Saxer2014_ancestorB %>%
  select(gene, population, frequency)

Saxer2014_ancestorB <- Saxer2014_ancestorB %>%
  replace_na(value = 0)

Saxer2014_ancestorB$frequency <- as.integer(as.logical(Saxer2014_ancestorB$frequency))

Saxer2014_ancestorB$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Saxer2014_ancestorB$gene)

write_csv(Saxer2014_ancestorB, here("data_in", "for_func", "Saxer2014_ancestorB.csv"))

ERConstraint(
  paper = "Saxer2014",
  dataset_name = "Saxer2014_ancestorB",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "BHI",
  generations = "765", 
  selective_pressure = "Rich media (BBL BHI)",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )


## Saxer2014_ancestorC:
Saxer2014_ancestorC <- read_csv(
  here("data_in", "original & usable", "Saxer2014", "Saxer2014_ancestorC_usable.csv")
  )

Saxer2014_ancestorC <- clean_names(Saxer2014, case = "snake")

colnames(Saxer2014_ancestorC) <- tolower(colnames(Saxer2014_ancestorC))

Saxer2014_ancestorC_out <- c(
  grep(out_patterns_column_gene, Saxer2014_ancestorC$gene),
  grep(out_patterns_column_details, Saxer2014_ancestorC$details)
  )

if (length(Saxer2014_ancestorC_out) > 0) {   
  Saxer2014_ancestorC <- Saxer2014[-Saxer2014_ancestorC_out,] 
} 

Saxer2014_ancestorC <- Saxer2014_ancestorC %>%
  select(gene, population, frequency)

Saxer2014_ancestorC <- Saxer2014_ancestorC %>%
  replace_na(value = 0)

Saxer2014_ancestorC$frequency <- as.integer(as.logical(Saxer2014_ancestorC$frequency))

Saxer2014_ancestorC$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Saxer2014_ancestorC$gene)

write_csv(Saxer2014_ancestorC, here("data_in", "for_func", "Saxer2014_ancestorC.csv"))

ERConstraint(
  paper = "Saxer2014",
  dataset_name = "Saxer2014_ancestorC",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "BHI",
  generations = "765", 
  selective_pressure = "Rich media (BBL BHI)",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Sherlock2013
Sherlock2013 <- read_csv(
  here("data_in", "original & usable", "Sherlock2013", "Sherlock2013_usable.csv")
  )

Sherlock2013 <- clean_names(Sherlock2013, case = "snake")

colnames(Sherlock2013) <- tolower(colnames(Sherlock2013))

Sherlock2013$details <- tolower(Sherlock2013$details)

names(Sherlock2013) <- gsub("gene", "Gene", names(Sherlock2013))

names(Sherlock2013) <- gsub("g", "", names(Sherlock2013))

names(Sherlock2013) <- gsub("Gene", "gene", names(Sherlock2013))
Sherlock2013_out <- c(
  grep(out_patterns_column_gene, Sherlock2013$gene),
  grep(out_patterns_column_details, Sherlock2013$details)
  )

if (length(Sherlock2013_out) > 0) {
  Sherlock2013 <- Sherlock2013[-Sherlock2013_out,]
}

Sherlock2013 <- Sherlock2013 %>%
  select(gene, population, "7", "70", "133", "196", "266", "322", "385", "448")

Sherlock2013 <- Sherlock2013 %>%
  replace_na(value = 0)

Sherlock2013$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sherlock2013$gene)


write_csv(Sherlock2013, here("data_in", "for_func", "Sherlock2013.csv"))

ERConstraint(
  paper = "Sherlock2013",
  dataset_name = "Sherlock2013",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "YPD",
  generations = c("7", "70", "133", "196", "266", "322", "385", "448"), 
  selective_pressure = "YPD",
  species = "Sac",
  who_analyzed = "MH",
  ploidy = "haploid"
 )

# Sherlock2019
Sherlock2019 <- read_csv(
  here("data_in", "original & usable", "Sherlock2019", "Sherlock2019_usable.csv")
  )

Sherlock2019 <- clean_names(Sherlock2019, case = "snake")

colnames(Sherlock2019) <- tolower(colnames(Sherlock2019))

names(Sherlock2019) <- gsub("x", "", names(Sherlock2019))

Sherlock2019_out <- c(
  grep(out_patterns_column_gene, Sherlock2019$gene),
  grep(out_patterns_column_details, Sherlock2019$details)
  )

if (length(Sherlock2019_out) > 0) {   
  Sherlock2019 <- Sherlock2019[-Sherlock2019_out,] 
} 

Sherlock2019 <- Sherlock2019 %>%
  select(gene, population, "50", "100", "150", "200", "250", "300", "350", "400", "450", "500")

Sherlock2019 <- Sherlock2019 %>%
  replace_na(value = 0)

Sherlock2019$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Sherlock2019$gene)

write_csv(Sherlock2019, here("data_in", "for_func", "Sherlock2019.csv"))

ERConstraint(
  paper = "Sherlock2019",
  dataset_name = "Sherlock2019",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "Davis minimal medium", 
  generations = c("50", "100", "150", "200", "250", "300", "350", "400", "450", "500"),
  selective_pressure = "Constant, glucose-limited environment", 
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid"
  )

# Suzuki2017
Suzuki2017 <- read_csv(
  here("data_in", "original & usable", "Suzuki2017", "Suzuki2017_usable.csv")
  )

Suzuki2017 <- clean_names(Suzuki2017, case = "snake")

colnames(Suzuki2017) <- tolower(colnames(Suzuki2017))

Suzuki2017_out <- c(
  grep(out_patterns_column_gene, Suzuki2017$gene),
  grep(out_patterns_column_details, Suzuki2017$details)
  )

if (length(Suzuki2017_out) > 0) {   
  Suzuki2017 <- Suzuki2017[-Suzuki2017_out,] 
} 

Suzuki2017 <- Suzuki2017 %>%
  select(gene, population, selective_pressure, frequency)

Suzuki2017 <- Suzuki2017 %>%
  replace_na(value = 0)

Suzuki2017$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Suzuki2017$gene)

write_csv(Suzuki2017, here("data_in", "for_func", "Suzuki2017.csv"))

ERConstraint(
  paper = "Suzuki2017",
  dataset_name = "Suzuki2017",
  timepoint_pressure_info = "multiple",
  structure = "long",
  environment = "M9",
  days = "33",
  selective_pressure = c("AMK & CP", "AMK & ENX", "CP & ENX", "AMK", "CP", "ENX"), 
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Tack2018
## Tack2018_control:
Tack2018_control <- read_csv(
  here("data_in", "original & usable", "Tack2018", "Tack2018_control_usable.csv")
  )

Tack2018_control <- clean_names(Tack2018_control, case = "snake")

colnames(Tack2018_control) <- tolower(colnames(Tack2018_control))

Tack2018_control_out <- c(
  grep(out_patterns_column_gene, Tack2018_control$gene),
  grep(out_patterns_column_details, Tack2018_control$details)
  )

if (length(Tack2018_control_out) > 0) {   
  Tack2018_control <- Tack2018_control[-Tack2018_control_out,] 
} 

Tack2018_control <- Tack2018_control %>%
  select(gene, population, selective_pressure, frequency)

Tack2018_control <- Tack2018_control %>%
  replace_na(value = 0)

Tack2018_control$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Tack2018_control$gene)

write_csv(Tack2018_control, here("data_in", "for_func", "Tack2018_control.csv"))

ERConstraint(
  paper = "Tack2018",
  dataset_name = "Tack2018_control",
  timepoint_pressure_info = "multiple",
  structure = "long",
  environment = "Rich defined media",
  generations = "2000",
  selective_pressure = c("RDM20", "RDM19", "RDM13"), 
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

## Tack2018_addicted:
Tack2018_addicted <- read_csv(
  here("data_in", "original & usable", "Tack2018", "Tack2018_addicted_usable.csv")
  )

Tack2018_addicted <- clean_names(Tack2018_addicted, case = "snake")

colnames(Tack2018_addicted) <- tolower(colnames(Tack2018_addicted))

Tack2018_addicted_out <- c(
  grep(out_patterns_column_gene, Tack2018_addicted$gene),
  grep(out_patterns_column_details, Tack2018_addicted$details)
  )

if (length(Tack2018_addicted_out) > 0) {   
  Tack2018_addicted <- Tack2018_addicted[-Tack2018_addicted_out,] 
} 

Tack2018_addicted <- Tack2018_addicted %>%
  select(gene, population, selective_pressure, frequency)

Tack2018_addicted <- Tack2018_addicted %>%
  replace_na(value = 0)

Tack2018_addicted$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Tack2018_addicted$gene)

write_csv(Tack2018_addicted, here("data_in", "for_func", "Tack2018_addicted.csv"))

ERConstraint(
  paper = "Tack2018",
  dataset_name = "Tack2018_addicted",
  timepoint_pressure_info = "multiple",
  structure = "long",
  environment = "Rich defined media",
  generations = "2000",
  selective_pressure = c("RDM20", "RDM19", "RDM13"), 
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Tenaillon2012
Tenaillon2012 <- read_csv(
  here("data_in", "original & usable", "Tenaillon2012", "Tenaillon2012_usable.csv")
  )

Tenaillon2012 <- clean_names(Tenaillon2012, case = "snake")

colnames(Tenaillon2012) <- tolower(colnames(Tenaillon2012))

Tenaillon2012_out <- c(
  grep(out_patterns_column_gene, Tenaillon2012$gene),
  grep(out_patterns_column_details, Tenaillon2012$details)
  )

if (length(Tenaillon2012_out) > 0) {   
  Tenaillon2012 <- Tenaillon2012[-Tenaillon2012_out,] 
} 

Tenaillon2012 <- gather(
  Tenaillon2012,
  population,
  frequency,
  "a0_f1_i1_r1": "a143_f1_i1_r1",
  factor_key = TRUE )

Tenaillon2012 <- Tenaillon2012 %>%
  select(gene, population, frequency)

Tenaillon2012 <- Tenaillon2012 %>%
  replace_na(value = 0)

Tenaillon2012$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Tenaillon2012$gene)

write_csv(Tenaillon2012, here("data_in", "for_func", "Tenaillon2012.csv"))

ERConstraint(
  paper = "Tenaillon2012",
  dataset_name = "Tenaillon2012",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Davis minimal medium",
  generations = "2000", 
  selective_pressure = "heat (42.2 degrees C)",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Tenaillon2016   
Tenaillon2016 <- read_csv(
  here("data_in", "original & usable", "Tenaillon2016", "Tenaillon2016_usable.csv")
  )
Tenaillon2016 <- clean_names(Tenaillon2016, case = "snake")

colnames(Tenaillon2016) <- tolower(colnames(Tenaillon2016))

Tenaillon2016_out <- c(
  grep(out_patterns_column_gene, Tenaillon2016$gene),
  grep(out_patterns_column_details, Tenaillon2016$details)
  )

if (length(Tenaillon2016_out) > 0) {
  Tenaillon2016 <- Tenaillon2016[-Tenaillon2016_out,]
}

Tenaillon2016 = Tenaillon2016 %>% tidyr::replace_na(
  list(
    x500_i1_r1 = 0,
    x500_i2_r1 = 0,
    x1000_i1_r1 = 0,
    x1000_i2_r1 = 0,
    x1500_i1_r1 = 0,
    x1500_i2_r1 = 0,
    x2000_i1_r1 = 0,
    x2000_i2_r1 = 0,
    x5000_i1_r1 = 0,
    x5000_i2_r1 = 0,
    x10000_i1_r1 = 0,
    x10000_i2_r1 = 0,
    x15000_i1_r1 = 0,
    x15000_i2_r1 = 0,
    x20000_i1_r1 = 0,
    x20000_i2_r1 = 0,
    x30000_i1_r1 = 0,
    x30000_i2_r1 = 0,
    x40000_i1_r1 = 0,
    x40000_i2_r1 = 0,
    x50000_i1_r1 = 0,
    x50000_i2_r1 = 0
  )
) 

Tenaillon2016 <- Tenaillon2016 %>% 
  transmute(gene,
            population,
            "500" = `x500_i1_r1`+`x500_i2_r1`,
            "1000" =`x1000_i1_r1`+`x1000_i2_r1`,
            "1500" =`x1500_i1_r1`+`x1500_i2_r1`,
            "2000"= `x2000_i1_r1`+`x2000_i2_r1`,
            "5000"=`x5000_i1_r1`+`x5000_i2_r1`,
            "10000"= `x10000_i1_r1`+`x10000_i2_r1`,
            "15000"=`x15000_i1_r1`+`x15000_i2_r1`,
            "20000"=`x20000_i1_r1`+`x20000_i2_r1`,
            "30000"= `x30000_i1_r1`+`x30000_i2_r1`,
            "40000"=`x40000_i1_r1`+`x40000_i2_r1`,
            "50000"=`x50000_i1_r1`+`x50000_i2_r1`
            ) 

Tenaillon2016$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Tenaillon2016$gene)

write_csv(Tenaillon2016, here("data_in", "for_func", "Tenaillon2016.csv"))

ERConstraint(
  paper = "Tenaillon2016",
  dataset_name = "Tenaillon2016",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "Davis minimal medium", 
  generations = c("500", "1000", "1500", "2000", '5000', '10000', '15000', '20000', '30000', '40000','50000'),
  selective_pressure = "Davis minimal medium", 
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid"
  )

# Wang2010
Wang2010 <- read_csv(
  here("data_in", "original & usable", "Wang2010", "Wang2010_usable.csv")
  )

Wang2010 <- clean_names(Wang2010, case = "snake")

colnames(Wang2010) <- tolower(colnames(Wang2010))

Wang2010_out <- c(
  grep(out_patterns_column_gene, Wang2010$gene),
  grep(out_patterns_column_details, Wang2010$details)
  )

if (length(Wang2010_out) > 0) {   
  Wang2010 <- Wang2010[-Wang2010_out,] 
} 
Wang2010 <- Wang2010 %>%
  select(gene, population, frequency)
Wang2010 <- Wang2010 %>%
  replace_na(value = 0)

Wang2010$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Wang2010$gene)

write_csv(Wang2010, here("data_in", "for_func", "Wang2010.csv"))

ERConstraint(
  paper = "Wang2010",
  dataset_name = "Wang2010",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "T-salts minimal medium",
  days = "37", 
  selective_pressure = "phosphate limitation",
  species = "Ecoli_K12",
  who_analyzed = "TL",
  ploidy = "haploid"
  )

# Wannier2018
Wannier2018 <- read_csv(
  here("data_in", "original & usable", "Wannier2018", "Wannier2018_usable.csv") 
  )

Wannier2018 <- clean_names(Wannier2018, case = "snake")

colnames(Wannier2018) <- tolower(colnames(Wannier2018))

Wannier2018_out <- c(
  grep(out_patterns_column_gene, Wannier2018$gene),
  grep(out_patterns_column_details, Wannier2018$details)
  )

if (length(Wannier2018_out) > 0) {   
  Wannier2018 <- Wannier2018[-Wannier2018_out,] 
} 

Wannier2018 = Wannier2018 %>% tidyr::replace_na(
  list(
    a1_f1_i1_r1 = 0,
    a1_f1_i2_r1 = 0,
    a2_f1_i1_r1 = 0,
    a2_f1_i2_r1 = 0,
    a3_f1_i1_r1 = 0,
    a3_f1_i2_r1 = 0, 
    a4_f1_i1_r1 = 0,
    a4_f1_i2_r1 = 0,
    a5_f1_i1_r1 = 0,
    a5_f1_i2_r1 = 0,
    a6_f1_i1_r1 = 0,
    a6_f1_i2_r1 = 0,  
    a7_f1_i1_r1 = 0,
    a7_f1_i2_r1 = 0,
    a8_f1_i1_r1 = 0,
    a8_f1_i2_r1 = 0,
    a9_f1_i1_r1 = 0,
    a9_f1_i2_r1 = 0, 
    a10_f1_i1_r1 = 0,
    a10_f1_i2_r1 = 0,
    a11_f1_i1_r1 = 0,
    a11_f1_i2_r1 = 0,
    a12_f1_i1_r1 = 0,
    a12_f1_i2_r1 = 0,  
    a13_f1_i1_r1 = 0,
    a13_f1_i2_r1 = 0,
    a14_f1_i1_r1 = 0,
    a14_f1_i2_r1 = 0
  )
) 

Wannier2018 <- Wannier2018 %>%
  transmute(
    gene = gene,
    a1 = `a1_f1_i1_r1`+ `a1_f1_i2_r1`,
    a2 = `a2_f1_i1_r1`+ `a2_f1_i2_r1`,
    a3 = `a3_f1_i1_r1`+ `a3_f1_i2_r1`, 
    a4 = `a4_f1_i1_r1`+ `a4_f1_i2_r1`,
    a5 = `a5_f1_i1_r1`+ `a5_f1_i2_r1`,
    a6 = `a6_f1_i1_r1`+ `a6_f1_i2_r1`,  
    a7 = `a7_f1_i1_r1`+ `a7_f1_i2_r1`,
    a8 = `a8_f1_i1_r1`+ `a8_f1_i2_r1`,
    a9 = `a9_f1_i1_r1`+ `a9_f1_i2_r1`, 
    a10 = `a10_f1_i1_r1`+ `a10_f1_i2_r1`,
    a11 = `a11_f1_i1_r1`+ `a11_f1_i2_r1`,
    a12 = `a12_f1_i1_r1`+ `a12_f1_i2_r1`,  
    a13 = `a13_f1_i1_r1`+ `a13_f1_i2_r1`,
    a14 = `a14_f1_i1_r1`+ `a14_f1_i2_r1`
    ) %>%
  gather(population, frequency, "a1" : "a14", factor_key = TRUE) %>%
  select(gene, population, frequency)

Wannier2018$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Wannier2018$gene)

write_csv(Wannier2018, here("data_in", "for_func", "Wannier2018.csv"))

ERConstraint(
  paper = "Wannier2018",
  dataset_name = "Wannier2018",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "M9 minimal medium",
  generations = "1100", 
  selective_pressure = "M9 minimal medium",
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid",
  strain_info = " All 321 UAG stop codons replaced with UAA; prfA deleted"
  )

# Wenger2011
Wenger2011 <- read_csv(
  here("data_in", "original & usable", "Wenger2011", "Wenger2011_usable.csv")
  )

Wenger2011 <- clean_names(Wenger2011, case = "snake")

colnames(Wenger2011) <- tolower(colnames(Wenger2011))

Wenger2011_out <- c(
  grep(out_patterns_column_gene, Wenger2011$gene),
  Wenger2011_out <-grep(out_patterns_column_details, Wenger2011$details)
  )

if (length(Wenger2011_out) > 0) {   
  Wenger2011 <- Wenger2011[-Wenger2011_out,] 
} 

Wenger2011 <- Wenger2011 %>%
  select(gene, population, frequency)

Wenger2011 <- Wenger2011 %>%
  replace_na(value = 0)

Wenger2011$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Wenger2011$gene)

write_csv(Wenger2011, here("data_in", "for_func", "Wenger2011.csv"))

ERConstraint(
  paper = "Wenger2011",
  dataset_name = "Wenger2011",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "YEPD (2% glucose)",
  generations = "250", 
  selective_pressure = "Continuous aerobic glucose limitation",
  species = "Sac",
  who_analyzed = "TL",
  ploidy = "diploid"
  )


# Wielgoss2016 
## Wielgoss2016_mucoid
Wielgoss2016_mucoid <- read_csv(
  here("data_in", "original & usable", "Wielgoss2016", "Wielgoss2016_mucoid_usable.csv")
  )

Wielgoss2016_mucoid <- clean_names(Wielgoss2016_mucoid, case = "snake")

colnames(Wielgoss2016_mucoid) <- tolower(colnames(Wielgoss2016_mucoid))

names(Wielgoss2016_mucoid) <- gsub("x", "", names(Wielgoss2016_mucoid))

Wielgoss2016_mucoid_out <- c(
  grep(out_patterns_column_gene, Wielgoss2016_mucoid$gene),
  grep(out_patterns_column_details, Wielgoss2016_mucoid$details)
  )

if (length(Wielgoss2016_mucoid_out) > 0) {   
  Wielgoss2016_mucoid <- Wielgoss2016_mucoid[-Wielgoss2016_mucoid_out,] 
}

Wielgoss2016_mucoid <- Wielgoss2016_mucoid %>%
  transmute(gene, population, frequency)

Wielgoss2016_mucoid <- Wielgoss2016_mucoid %>%
  replace_na(value = 0)

Wielgoss2016_mucoid$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Wielgoss2016_mucoid$gene)

write_csv(Wielgoss2016_mucoid, here("data_in", "for_func", "Wielgoss2016_mucoid.csv"))

ERConstraint(
  paper = "Wielgoss2016",
  dataset_name = "Wielgoss2016_mucoid",
  timepoint_pressure_info = "single",
  structure = "long",
  environment = "Lysogeny broth",
  generations = "10", 
  selective_pressure = "Lysogeny broth",
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid",
  strain_info = "mucoid"
  )

## Wielgoss2016 non-Mucoid
Wielgoss2016_nonMucoid <- read.csv(
  here("data_in", "original & usable", "Wielgoss2016", "Wielgoss2016_nonMucoid_usable.csv"),
  encoding = "UTF-8"
  )

Wielgoss2016_nonMucoid <- clean_names(Wielgoss2016_nonMucoid, case = "snake")

colnames(Wielgoss2016_nonMucoid) <- tolower(colnames(Wielgoss2016_nonMucoid))

names(Wielgoss2016_nonMucoid) <- gsub("x", "", names(Wielgoss2016_nonMucoid))

Wielgoss2016_nonMucoid_out <- c(
  grep(out_patterns_column_gene, Wielgoss2016_nonMucoid$gene),
  grep(out_patterns_column_details, Wielgoss2016_nonMucoid$details)
  )

if (length(Wielgoss2016_nonMucoid_out) > 0) {   
  Wielgoss2016_nonMucoid <- Wielgoss2016_nonMucoid[-Wielgoss2016_nonMucoid_out,] 
} 

Wielgoss2016_nonMucoid <- Wielgoss2016_nonMucoid %>%
  select(gene, population, "2", "10")

Wielgoss2016_nonMucoid <- Wielgoss2016_nonMucoid %>%
  replace_na(value = 0)

Wielgoss2016_nonMucoid$gene <- gsub("[^[:alnum:][:blank:]&/\\-]", "", Wielgoss2016_nonMucoid$gene)

write_csv(Wielgoss2016_nonMucoid, here("data_in", "for_func", "Wielgoss2016_nonMucoid.csv"))

ERConstraint(
  paper = "Wielgoss2016",
  dataset_name = "Wielgoss2016_nonMucoid",
  timepoint_pressure_info = "multiple",
  structure = "wide",
  environment = "Lysogeny broth",
  generations = c("2", "10"), 
  selective_pressure = "Lysogeny broth",
  species = "Ecoli_K12",
  who_analyzed = "MH",
  ploidy = "haploid",
  strain_info = "non-mucoid"
  )
