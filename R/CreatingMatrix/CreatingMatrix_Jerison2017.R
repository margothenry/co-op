#####################
#Run every time
library(tidyverse)
library(readr)
library(devtools)
library(dgconstraint)
geneNumbers <- read_csv("/Users/Margot/Desktop/Co-op/data-in/GeneDatabase.csv")

#################################################################################################
#Read in the data and type in specifics
#################################################################################################
data <- read_csv("/Users/Margot/Desktop/Co-op/data-in/Jerison 2017.csv")
Population <- data$Population
paper <- "Jersion2017"
environment <- "Synthetic-HT-high temperature"
selective_pressure <- unique(data$selective_pressure)

#"Sac" or "Ecoli_K12" or "Ecoli_O157-H7" 
numGenes <- filter(geneNumbers, Species == "Sac")$NumGenes
#generations <- 
#TRUE OF FALSE
multiple <- TRUE
matrix <-FALSE
medium <-  TRUE
#if need to filter out mutations(Not SNPs) or anything else(like intergenics) from data do so here
data <-data %>% filter(Gene != "Intergenic") %>% filter(selective_pressure == "OT")
Founder <- unique(data$Founder)
###############################################################################################

if (multiple == TRUE && matrix == FALSE && medium == FALSE) {
  numLineages <- c()
  num_parallel_genes <- c()
  num_non_parallel_genes <- c()
  parallel_genes <- c()
  c_hyper <- c()
  p_chisq <- c()
  estimate <- c()
  
  for (g in generations) {
   #####################################################################################################################
     print(g)
    data.1 <- data %>% 
      arrange(Gene) %>%
      drop_na(Gene) %>%
      drop_na(Population)
   
    data.g <- data.1 %>% 
      select(Population, Gene, Prop = g) %>% 
      filter(Prop >0)
    ##################################################################################################################
    
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
      num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), Genes = row.names(data.array))
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
    parallel_genes <- append(parallel_genes, paste0(genes_parallel$Genes, collapse=", ")) 
    
    
    full_matrix <- rbind(data.array, array(0,c(numGenes-total_genes,ncol(data.array))))
    
    #write_csv(data.g, path = paste0("data-out/paper_", g, ".csv"))
    
    c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
    p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
    estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))
  }
  
} else if(multiple == FALSE && matrix == FALSE && medium == FALSE) {
  data.1 <- data %>%
    arrange(Gene) %>%
    drop_na(Gene) %>%
    drop_na(Population)
  
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
  
  
} else if(multiple == TRUE && medium == TRUE && matrix == FALSE){
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
    drop_na(selective_pressure) %>%
    select(selective_pressure,Founder, Gene, frequency, Clone) 
  
  for(j in Founder) {
    print(j)
    data.j <- data.1 %>% 
    filter(Founder == j)
    
    num_genes <- length((unique(data.j$Gene)))
    num_lineages <- length(unique(data.j$Clone))
    data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.j$Gene), unique(data.j$Clone)))
    
    for(i in 1:num_lineages) {
      sub <- subset(data.j, data.j$Clone == unique(data.j$Clone)[i])
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
}else{
  data.1 <- data %>%
    arrange(Gene) %>%
    drop_na(Gene) %>%
    drop_na(Population)
  data.2 <- data.1 %>%
    select(Population)

  data.matrix <- as.matrix(data.2)
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

}

#######if NOT MULTIPLE#######
df <- tibble( paper = paper, environment, temp = temp, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)

#######if MULTIPLE for generation#######
df <- tibble( paper = paper, environment = environment, generations, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)

#######if MULTIPLE for medium(Temperature, environment...)#######
df <- tibble( paper = paper, environment = environment,Founder, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)

#######SAVE ANALYSIS######
write_csv(df, path = paste0("/Users/Margot/Desktop/Co-op/data-out/", paper, "_Analysis_Temp.csv"))
