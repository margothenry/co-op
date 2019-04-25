#####################
#Run every time
#####################
library(tidyverse)
library(readr)
library(readr)
library(devtools)
library(dgconstraint)
geneNumbers <- read_csv("data-in/GeneDatabase.csv")

###############################################################
#Read in the data and (for now) type in specifics
###############################################################
data <- read_csv("data-in/Sandberg_2014.csv")
population <- c("ALE1", "ALE2", "ALE3", "ALE4", "ALE5", "ALE6", "ALE7", "ALE8", "ALE9", "ALE10")
paper <- "Sandberg2014"
environment <- "glucose minimal media"
numGenes <- filter(geneNumbers, Species == "Ecoli_K12")$NumGenes
multipe <- FALSE

###############################################################
#none of this is specific! Run.
###############################################################

###############################################################
#Organize the data in a tidy manner
###############################################################
data.1 <- data %>% 
  arrange(Gene) %>%
  drop_na(Gene) %>%
  drop_na(population) 
data.2 <- data.1 %>%
  select(population)

data.matrix <- as.matrix(data.2)
rownames(data.matrix) <-data.1$Gene

num_parallel <- data.frame(data.matrix, Count=rowSums(data.matrix, na.rm = FALSE, dims = 1), Gene = row.names(data.matrix))

genes_parallel <- num_parallel %>% 
  as_tibble() %>% 
  filter(Count > 1)

Non_genes_parallel <- num_parallel %>% 
  as_tibble() %>% 
  filter(Count == 1)

###############################################################
#Calculation of numbers of genes of different classes & stats
###############################################################
if(multiple = TRUE){
  num_parallel_genes <- c()
  num_non_parallel_genes <- c()
  parallel_genes <- c()
  c_hyper <- c()
  p_chisq <- c()
  estimate <- c()
  num_parallel_genes <- append(num_parallel_genes, nrow(genes_parallel))
  num_non_parallel_genes <- append(num_non_parallel_genes, nrow(Non_genes_parallel))
  parallel_genes <- append(parallel_genes, paste0(genes_parallel$Gene, collapse=", "))
  full_matrix <- rbind(data.matrix, array(0,c(numGenes-total_genes,ncol(data.matrix))))
  
  c_hyper <- append(c_hyper, pairwise_c_hyper(data.matrix))
  p_chisq <- append(p_chisq, allwise_p_chisq(data.matrix, num_permute = 200))
  estimate <- append(estimate, estimate_pa(data.matrix,ndigits = 4, show.plot = T))
} else{
  num_parallel_genes <- nrow(genes_parallel)
  num_non_parallel_genes <- nrow(Non_genes_parallel)
  parallel_genes <- paste0(genes_parallel$Gene, collapse=", ")
  total_genes <- num_non_parallel_genes + num_parallel_genes
  
  full_matrix <- rbind(data.matrix, array(0,c(numGenes-total_genes,ncol(data.matrix)))) 
  
  c_hyper <- pairwise_c_hyper(full_matrix)
  p_chisq <- allwise_p_chisq(full_matrix, num_permute = 200)
  estimate <- estimate_pa(full_matrix, ndigits = 4, show.plot = T)
}

df <- tibble(paper = paper, environment= environment, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
write_csv(df, path = paste0("data-out/", paper, "_Analysis.csv"))
