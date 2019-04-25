#####################
#Run every time
library(tidyverse)
library(readr)
library(devtools)
library(dgconstraint)
geneNumbers <- read_csv("/Users/Margot/Desktop/Co-op/data-in/GeneDatabase.csv")

#################################################################################################
#Read in the data and (for now) type in specifics
#################################################################################################
data <- read_csv("/Users/Margot/Desktop/Co-op/data-in/Sandberg2016.csv")
population <- c("ALE1", "ALE2", "ALE3", "ALE4", "ALE5", "ALE6")
paper <- "Sandberg2016"
environment <- " Davids minimal medium"
#"Sac" or "Ecoli_K12" or "Ecoli_O157-H7" 
numGenes <- filter(geneNumbers, Species == "Ecoli_K12")$NumGenes
generations <-c("Flask 23", "Flask 58", "Flask 133") 
multiple <- TRUE

if (multiple == TRUE) {
  numLineages <- c()
  num_parallel_genes <- c()
  num_non_parallel_genes <- c()
  parallel_genes <- c()
  c_hyper <- c()
  p_chisq <- c()
  estimate <- c()
  
  for (g in generations) {
   #####################################################################################################################
    #next block might change depending on names on data, which data you need to filter out, generations or environments
     print(g)
    data.1 <- data %>% 
      arrange(Gene) %>%
      drop_na(Gene) %>%
      drop_na(Population)
   
    data.g <- data.1 %>% 
      select(Population, Gene, Prop = g) %>% 
      filter(Prop >0)
    ####################################################################################################################
    
    #number of unique genes
    num_genes <- length((unique(data.g$Gene)))
    
    #number of lineages
    num_lineages <- length(unique(data.g$Population))
    numLineages <- append(numLineages, num_lineages)
    
    #making the matrix
    data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.g$Gene), unique(data.g$Population)))
    
    for(i in 1:num_lineages) {
      sub <- subset(data.g, data.g$Population == unique(data.g$Population)[i])
      sub2 <- subset(sub, Prop > 0)
      geneRows <- which(row.names(data.array) %in% sub2$Gene)
      data.array[geneRows, i] <- 1
      num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), Gene = row.names(data.array))
    }
    
    genes_parallel <- num_parallel %>% 
      as_tibble() %>% 
      filter(Count > 1)
    num_parallel_genes_g <- nrow(genes_parallel)
    
    Non_genes_parallel <- num_parallel %>% 
      as_tibble() %>% 
      filter(Count == 1)
    num_non_parallel_genes_g <- nrow(Non_genes_parallel)
    total_genes <- num_non_parallel_genes_g + num_parallel_genes_g
    
    num_parallel_genes <- append(num_parallel_genes, num_parallel_genes_g)
    num_non_parallel_genes <- append(num_non_parallel_genes, num_non_parallel_genes_g)
    parallel_genes <- append(parallel_genes, paste0(genes_parallel$Gene, collapse=", ")) 
    
    
    full_matrix <- rbind(data.array, array(0,c(numGenes-total_genes,ncol(data.array))))
    
    write_csv(data.g, path = paste0("data-out/Sandberg2016_", g, ".csv"))
    
    c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
    p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
    estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))
  }
  
} else {
  #####################################################################################################################
  #next block might change depending on names on data, which data you need to filter out, generations or environments
  data.1 <- data %>%
    arrange(Gene) %>%
    drop_na(Gene) %>%
    drop_na(population)
  data.2 <- data.1 %>%
    select(population)
  #####################################################################################################################

  data.matrix <- as.matrix(data.2)
  rownames(data.matrix) <-data.1$Gene

  num_parallel <- data.frame(data.matrix, Count=rowSums(data.matrix, na.rm = FALSE, dims = 1), gene = row.names(data.matrix))

  genes_parallel <- num_parallel %>%
    as_tibble() %>%
    filter(Count > 1)


  Non_genes_parallel <- num_parallel %>%
    as_tibble() %>%
    filter(Count == 1)


  num_parallel_genes <- nrow(genes_parallel)
  num_non_parallel_genes <- nrow(Non_genes_parallel)
  total_genes <- num_non_parallel_genes + num_parallel_genes
  parallel_genes <- paste0(genes_parallel$gene, collapse=", ")

  full_matrix <- rbind(data.matrix, array(0,c(numGenes-total_genes,ncol(data.matrix))))


  c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
  p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
  estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))

}

#######if FALSE#######
df <- tibble( paper = paper, environment = environment, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
write_csv(df, path = paste0("/Users/Margot/Desktop/Co-op/data-out/", paper, "_Analysis.csv"))
#######if TRUE#######
df <- tibble( paper = paper, environment = environment, generations, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)