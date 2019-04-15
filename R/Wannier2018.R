#' Calculations for a Single Generation Functionin Matrix Form
#'
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a data set with a single generation.
#' 
#' @param paper the data in csv that you want to analyze
#' @param environment The environment in which the experiment occured
#' @param species Specify if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7" 
#' @param population a list of populations in the data
#' @return a table with all the calculated infromation
#' @export 
#' @examples 
#' singlematrix("Author2018","YPD", "Sac", c("P1", "P2", "P3" ,"P4", "P5"))
#'
#'
#####################

singlematrix <- function(paper, environment, species, population)
  
  #Run every time - read in descrption file for package
library(tidyverse)
library(readr)
library(devtools)
library(dgconstraint)
library(Hmisc)

geneNumbers <- read_csv(file.path(getwd(),"data-in/GeneDatabase.csv"))
#########################################################################################
#Read in the data and type in specifics
#########################################################################################
paper <- "Wannier2018"
data <- read_csv(file.path(getwd(), "data-in", paste0(paper, ".csv")))
environment <- "Glucose Mininal Medium"
data <-data %>% filter(Details != "intergenic") 
data<- data %>% mutate(A1= `A1 F1 I1 R1`+ `A1 F1 I2 R1`, A2 = `A2 F1 I1 R1`+ `A2 F1 I2 R1`, A3 = `A3 F1 I1 R1`+ `A3 F1 I2 R1`, A4 = `A4 F1 I1 R1`+ `A4 F1 I2 R1`, A5 = `A5 F1 I1 R1`+ `A5 F1 I2 R1`, A6 = `A6 F1 I1 R1`+ `A6 F1 I2 R1`, A7 = `A7 F1 I1 R1`+ `A7 F1 I2 R1`, A8 = `A8 F1 I1 R1`+ `A8 F1 I2 R1`, A9 = `A9 F1 I1 R1`+ `A9 F1 I2 R1`, A10 = `A10 F1 I1 R1`+ `A10 F1 I2 R1`,A11 = `A1 F1 I1 R1`+ `A11 F1 I2 R1`, A12 = `A12 F1 I1 R1`+ `A12 F1 I2 R1`, A13 = `A13 F1 I1 R1`+ `A13 F1 I2 R1`, A14 = `A14 F1 I1 R1`+ `A14 F1 I2 R1`) %>% replace(is.na(.), 0) 
population <- c("A1" , "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14")
  #have users manually input species number if not using these species
  #if(species %nin% c("Sac", "Ecoli_K12", "Ecoli_O157-H7" )){
species <-"Ecoli_K12"
numGenes <- filter(geneNumbers, Species == species)$NumGenes


{
  numLineages <- c()
  num_parallel_genes <- c()
  num_non_parallel_genes <- c()
  parallel_genes <- c()
  c_hyper <- c()
  p_chisq <- c()
  estimate <- c()
  
  data.1 <- data %>%
  arrange(Gene) %>%
  drop_na(Gene) 
data.2 <- data.1 %>%
  select(population)

data.matrix <- as.matrix(data.2)
data.matrix[data.matrix > 0]<-1

rownames(data.matrix) <-data.1$Gene

num_parallel <- data.frame(data.matrix, Count=rowSums(data.matrix, na.rm = FALSE, dims = 1), Genes = row.names(data.matrix))

genes_parallel <- num_parallel %>%
  as_tibble() %>%
  filter(Count > 1)


Non_genes_parallel <- num_parallel %>%
  as_tibble() %>%
  filter(Count == 1)


num_parallel_genes <- nrow(genes_parallel)
num_non_parallel_genes <- nrow(Non_genes_parallel)
total_genes <- num_non_parallel_genes + num_parallel_genes
parallel_genes <- paste0(genes_parallel$Genes, collapse=", ")

full_matrix <- rbind(data.matrix, array(0,c(numGenes-total_genes,ncol(data.matrix))))


c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))

c_hyper[c_hyper <= 0] <- 0
c_hyper[c_hyper == "NaN"] <- 0

df <- tibble( paper = paper, environment = environment, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)

newdir <- file.path(getwd(), "data-out")
if (!file.exists(newdir)){
  dir.create(newdir, showWarnings = FALSE)
  cat(paste("\n\tCreating new directory: ", newdir), sep="")
}

filename <- file.path(getwd(), "data-out", paste(paper, "_Analysis.csv", sep=""))
write.csv(df, file=filename, row.names=FALSE)

}
