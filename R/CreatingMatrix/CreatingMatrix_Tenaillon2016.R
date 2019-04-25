
library(tidyverse)
library(readr)
library(devtools)
library(dgconstraint)

paper <-"Tenaillon2016_original"
environment <-"David minimal medium"
species <-"Ecoli_K12"
generations <-c("g500", 'g1000', 'g1500', 'g2000', 'g5000', 'g10000', 'g20000', 'g30000', 'g40000','g50000')   


multigen_c_hyper <- function(paper, environment, species, numGenes = NA, generations){

geneNumbers <- read_csv(file.path(getwd(),"data-in/GeneDatabase.csv"))
data<- read_csv(file.path(getwd(),"data-in/Tenaillon2016_original.csv"))

data <- data %>% 
  transmute(Gene = data$Gene, Population = data$Population, Details = data$Details, `g500` = `500 I1 R1`+`500 I2 R1`, `g1000` = `1000 I1 R1`+`1000 I2 R1`, `g1500` =`1500 I1 R1`+`1500 I2 R1`, `g2000` = `2000 I1 R1`+`2000 I2 R1`, `g5000` = `5000 I1 R1`+`5000 I2 R1`, `g10000` = `10000 I1 R1`+`10000 I2 R1`, `g15000` = `15000 I1 R1`+`15000 I2 R1`, `g20000` = `20000 I1 R1`+`20000 I2 R1`, `g30000` = `30000 I1 R1`+`30000 I2 R1`, `g40000` = `40000 I1 R1`+`40000 I2 R1`, `g50000` = `50000 I1 R1`+`50000 I2 R1`) %>% replace(is.na(.), 0)

data<- data %>% 
  filter(Details != "intergenic")

if (species %in% geneNumbers$Species){
  numGenes <- filter(geneNumbers, Species == species)$NumGenes  
}

if(is.na(numGenes)){
  prompt <- "Your species is unspecified or not in our database. How many genes does it have? \n"
  numGenes <- as.numeric(readline(prompt))
}
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
  
  df <- tibble( paper = paper, environment = environment, c_hyper = round(c_hyper, 3),generations, p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
  
  filename2 <- file.path(getwd(), "data-out", paste(paper, "_Analysis.csv", sep=""))
  
  write.csv(df, file=filename2, row.names=FALSE)
  
}
