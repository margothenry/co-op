#Run every time
library(tidyverse)
library(readr)
library(devtools)
library(dgconstraint)
library(Hmisc)

geneNumbers <- read_csv(file.path(getwd(),"data-in/GeneDatabase.csv"))
#########################################################################################
#Read in the data and type in specifics
#########################################################################################
paper <- "Payen2016"
data <- read_csv(file.path(getwd(), "data-in", paste0(paper, ".csv")))
selective_pressure <- c("phosphate", "sulfate", "glucose")
#"Sac" or "Ecoli_K12" or "Ecoli_O157-H7" 
species <- "Sac"
numGenes <- filter(geneNumbers, Species == species)$NumGenes
#have users manually input species number if not using these species
#if(species %nin% c("Sac", "Ecoli_K12", "Ecoli_O157-H7" )){
#}

{
  data.1 <- data %>%
    arrange(Gene) %>%
    drop_na(Gene) %>%
    drop_na(Population) %>% 
    filter(frequency != "clone") %>% 
    drop_na(frequency)
  
  num_genes <- length((unique(data.1$Gene)))
  num_lineages <- length(unique(data.1$Population))
  
  data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.1$Gene), unique(data.1$Population)))
  
  for(i in 1:num_lineages) {
    sub <- subset(data.1, data.1$Population == unique(data.1$Population)[i])
    sub2 <- subset(sub, frequency > 0)
    geneRows <- which(row.names(data.array) %in% sub2$Gene)
    data.array[geneRows, i] <- 1
    num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), Genes = row.names(data.array))
  }
  
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
  
  full_matrix <- rbind(data.array, array(0,c(numGenes-total_genes,ncol(data.array))))
  
  
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
  
  filename <- file.path(getwd(), "data-out", paste(paper, "_multiple_Analysis.csv", sep=""))
  write.csv(df, file=filename, row.names=FALSE)
}
 

#multiple selective pressure
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
  drop_na(Gene) %>%
  drop_na(Population) %>% 
  filter(frequency != "clone") %>% 
  drop_na(frequency)

for(j in selective_pressure) {
  print(j)
  data.j <- data.1 %>% 
    filter(selective_pressure == j)
  
  num_genes <- length((unique(data.j$Gene)))
  num_lineages <- length(unique(data.j$Population))
  data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.j$Gene), unique(data.j$Population)))
  
  for(i in 1:num_lineages) {
    sub <- subset(data.j, data.j$Population == unique(data.j$Population)[i])
    sub2 <- subset(sub, frequency > 0)
    geneRows <- which(row.names(data.array) %in% sub2$Gene)
    data.array[geneRows, i] <- 1
    num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), Genes = row.names(data.array))
  }
  
  genes_parallel <- num_parallel %>% 
    as_tibble() %>% 
    filter(Count > 1)
  num_parallel_genes_j <- nrow(genes_parallel)
  
  Non_genes_parallel <- num_parallel %>% 
    as_tibble() %>% 
    filter(Count == 1)
  num_non_parallel_genes_j <- nrow(Non_genes_parallel)
  total_genes <- num_non_parallel_genes_j + num_parallel_genes_j
  
  num_parallel_genes <- append(num_parallel_genes, num_parallel_genes_j)
  num_non_parallel_genes <- append(num_non_parallel_genes, num_non_parallel_genes_j)
  parallel_genes <- append(parallel_genes, paste0(genes_parallel$Genes, collapse=", ")) 
  
  
  full_matrix <- rbind(data.array, array(0,c(numGenes-total_genes,ncol(data.array))))
  
  #write_csv(data.g, path = paste0("data-out/paper_", j, ".csv"))
  
  c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
  p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
  estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))
  c_hyper[c_hyper <= 0] <- 0
  c_hyper[c_hyper == "NaN"] <- 0
} 
  df <- tibble( paper = paper, selective_pressure, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
  
}
