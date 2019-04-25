#####################
#Run every time
library(tidyverse)
library(readr)
library(devtools)
library(dgconstraint)
library(Hmisc)

#################################################################################################
#Read in the data and (for now) type in specifics
#################################################################################################
data <- read_csv("/Users/Margot/Desktop/Co-op/data-in/original/Creamer2016_OG.csv")
population <- c("K0001","K0011","K0002","K0003","K0014","K0015","KB026","K0022","K0023","K0006","K0007","K0010", "K0020","K0019","K0031","K0030")
paper <- "Creamer2016"
environment <- "LBK medium"
species <- "Ecoli_K12"
geneNumbers <- read_csv(file.path(getwd(),"data-in/GeneDatabase.csv"))
collapseMutations = TRUE
data <- data %>% 
  filter(Annotation!= "intergenic") %>% 
  filter(!str_detect(Annotation,"insH")) %>% replace(is.na(.), 0)

if (species %in% geneNumbers$Species){
  numGenes <- filter(geneNumbers, Species == species)$NumGenes  
}

if(is.na(numGenes)){
  prompt <- "Your species is unspecified or not in our database. How many genes does it have? \n"
  numGenes <- as.numeric(readline(prompt))
}

data.1 <- data %>%
  arrange(Gene) %>%
  drop_na(Gene) %>%
  drop_na(population)

num_lineages <- length(unique(population))
num_genes <- length((unique(data.1$Gene)))
data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.1$Gene), unique(population)))

if(collapseMutations){
  multiple_entry_genes <- subset(table(data.1$Gene), table(data.1$Gene) >1)
  
  # These are our genes with ony a single mutation
  single_mutation_genes <- subset(data.1, Gene %nin% names(multiple_entry_genes))  
  single_mutation_genes <- single_mutation_genes[, c(unique(population), "Gene")]
  
  # These are the genes with multiple mutations. It may be the case in the future or in some circumstances that you want to know parallelism at the mutation rather than gene level. In that case don't include this.
  multiple_mutation_genes <- subset(data.1, Gene %in% names(multiple_entry_genes))  
  multi_genes_matrix <- data.frame(Gene = names(multiple_entry_genes))
  for (k in 1:length(multiple_entry_genes)){
    sub <- subset(data.1, Gene == names(multiple_entry_genes)[k])
    for (j in unique(population)){
      multi_genes_matrix[k, j] <- sum(sub[1:nrow(sub), j])
    }
  }
  
  data.1 <- rbind(single_mutation_genes, multi_genes_matrix)
  data.1 <- data.1 %>% 
    arrange(Gene) 
}


num_parallel <- data.frame(data.1[, population], Count=rowSums(data.1[, population], na.rm = FALSE, dims = 1),row.names= data.1$Gene)

genes_parallel <- num_parallel %>%
  as_tibble() %>%
  filter(Count > 1)

Non_genes_parallel <- num_parallel %>%
  as_tibble() %>%
  filter(Count == 1)

num_parallel_genes <- nrow(genes_parallel)
num_non_parallel_genes <- nrow(Non_genes_parallel)
total_genes <- num_non_parallel_genes + num_parallel_genes

extraGenes <- data.frame(array(0, dim = c(numGenes-total_genes, length(population))))
names(extraGenes) <- names(num_parallel[, population])
full_matrix <- rbind(num_parallel[, population], extraGenes)

c_hyper <- c()
p_chisq <- c()
estimate <- c()
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
