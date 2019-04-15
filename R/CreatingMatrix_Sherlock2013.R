library(tidyverse)
library(devtools)
library(dgconstraint)
library(ggplot2)

data <- read_csv("data-in/Sherlock2013.csv")
geneNumbers <- read_csv("data-in/GeneDatabase.csv")
paper <- "Sherlock2013"
environment <- "YPD"
numGenes <- filter(geneNumbers, Species == "Sac")$NumGenes


#define empty vectors
numLineages <- c()
num_parallel_genes <- c()
num_non_parallel_genes <- c()
parallel_genes <- c()
c_hyper <- c()
p_chisq <- c()
estimate <- c()

generations <- c("7","70", "133","196","266", "322","385","448")
names(data)[11:18] <- generations

 

#Loop over all available generations
for (g in generations){
  print(g)
  data.1 <- data %>% 
    arrange(Gene) %>%
    drop_na(Gene) %>%
    drop_na(Population)
  #create a new dataframe that filters the genes with a final proportion > 0 and only keep the columns we're intersted in for now
  data.g <- data.1 %>% 
    select(Population, Gene, Prop = g) %>% 
    filter(Prop >0)
  
  #number of unique genes
  num_genes <- length((unique(data.1$Gene)))
  
  #number of lineages
  num_lineages <- length(unique(data.g$Population))
  numLineages <- append(numLineages, num_lineages)
  #set up matrix
  #add labels  onto the matrix ("Lineage" are columns, "Gene" are the row)
  data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.1$Gene), unique(data.g$Population)))
  
  
  #for each lineage, update the column adding 1's into the rows that correspond to genes that were identified in the dataset
  #which row corresponds to the genes in sub
  #turn each into a 1
  for(i in 1:num_lineages) {
    sub <- subset(data.g, data.g$Population == unique(data.g$Population)[i])
    sub2 <- subset(sub, Prop > 0)
    geneRows <- which(row.names(data.array) %in% sub2$Gene)
    data.array[geneRows, i] <- 1
    num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), Gene = row.names(data.array))
  }
  
  #pull out from data how many genes are in > 1 lineage (row sums > 1)
  #pull out what those genes are
  
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
  
  write_csv(data.g, path = paste0("data-out/Sherlock2013_", g, ".csv"))
  
  c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
  p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
  estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))
}


df <- tibble( paper = paper, environment = environment, generations, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
write_csv(df, path = paste0("data-out/", paper, "_Analysis.csv"))
