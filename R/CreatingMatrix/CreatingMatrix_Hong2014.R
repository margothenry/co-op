#####################
#Run every time
#####################
library(tidyverse)
library(devtools)
library(dgconstraint)
geneNumbers <- read_csv("data-in/GeneDatabase.csv")

##############################################################
#Read in the data and (for now) type in specifics
##############################################################
data <- read_csv("data-in/Hong2014.csv")
population <- data$population
paper <- "Hong2014"
environment <- c("Ammonium", "Arginine", "Urea", "Allantoin")
#might not need this for this data set.. we'll see
numGenes <- filter(geneNumbers, Species == "Sac")$NumGenes

numLineages <- c()
num_parallel_genes <- c()
num_non_parallel_genes <- c()
parallel_genes <- c()
c_hyper <- c()
p_chisq <- c()
estimate <- c()

data.1 <- data %>% 
  filter( type != 'DIP', type !='WCA', type != 'IND', type != 'AMP', type != 'DEL') %>%
  arrange(Gene) %>%
  drop_na(Gene) %>%
  drop_na(population) %>%
  select(`selective pressure`,population, Gene, frequency)


for(j in environment) {
#filter so that it will run for each environment !!!rows<-select, columns<-select!!!
  print(j)
  data.j <- data.1 %>% 
    filter(`selective pressure` == j)  
 #change num_lineages and data.array so that it will change for each environment  
  num_genes <- length((unique(data.j$Gene)))
  num_lineages <- length(unique(data.j$population))
  data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.j$Gene), unique(data.j$population)))
  
  #chnage the for loop so that the data will change with the envrironment 
  for(i in 1:num_lineages) {
    sub <- subset(data.j, data.j$population == unique(data.j$population)[i])
    sub2 <- subset(sub, freq > 0)
    geneRows <- which(row.names(data.array) %in% sub2$Gene)
    data.array[geneRows, i] <- 1
    num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), Genes = row.names(data.array))
  }
  #same as before
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
}

df <- tibble( paper = paper, environment, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
write_csv(df, path = paste0("/Users/Margot/Desktop/Co-op/data-out/", paper, "_Analysis.csv"))
