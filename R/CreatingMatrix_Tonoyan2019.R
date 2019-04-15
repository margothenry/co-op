
#Run every time - read in descrption file for package
library(tidyverse)
library(readr)
library(devtools)
library(dgconstraint)
library(Hmisc)

geneNumbers <- read_csv(file.path(getwd(),"data-in/GeneDatabase.csv"))
########################################################################################
#Read in the data and type in specifics
########################################################################################
paper <- "Tonoyan_2019"
species <- "Ecoli_ATCC"
numGenes <- filter(geneNumbers, Species == "Ecoli_ATCC")$NumGenes
data <- read_csv(file.path(getwd(), "data-in", paste0(paper, ".csv")))
environment <- "LB medium"
generations <- c("1", "14", "20")

data <-data %>% filter(Gene != "Intergenic") 

##################################################################################################

{ 
  numLineages <- c()
  num_parallel_genes <- c()
  num_non_parallel_genes <- c()
  parallel_genes <- c()
  c_hyper <- c()
  p_chisq <- c()
  estimate <- c()
  
  for (g in generations) {
    
    print(g)
    data.1 <- data %>% 
      arrange(Gene) %>%
      drop_na(Gene) %>%
      drop_na(Population)
    
    data.g <- data.1 %>% 
      select(Population, Gene, Prop = g) %>% 
      filter(Prop >0)
    
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
    
    newdir <- file.path(getwd(), "data-out")
    if (!file.exists(newdir)){
      dir.create(newdir, showWarnings = FALSE)
      cat(paste("\n\tCreating new directory: ", newdir), sep="")
    }
    
    filename1 <- file.path(getwd(), "data-out", paste0("/", paper, "_", g, ".csv"))
    write.csv(data.g, file=filename1, row.names=FALSE)
    
    c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
    p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
    estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))
    
    c_hyper[c_hyper <= 0] <- 0
    c_hyper[c_hyper == "NaN"] <- 0
  }
  
  df <- tibble( paper = paper, environment = environment, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
  
  filename2 <- file.path(getwd(), "data-out", paste(paper, "_Analysis.csv", sep=""))
  
  write.csv(df, file=filename2, row.names=FALSE)
  
}
