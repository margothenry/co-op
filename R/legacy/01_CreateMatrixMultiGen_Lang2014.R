library(tidyverse)
#Installing dgconstraint package
library(devtools)
library(dgconstraint)
library(ggplot2)
geneNumbers <- read_csv("~/Desktop/Co-op/data-in/GeneDatabase.csv")
#Read in the data
Lang2014 <- read_csv("~/Desktop/Co-op/data-in/Lang2014Nature.csv")
generations <- c("0","140","240","335",'415','505','585','665','745','825','910','1000')
paper <- "Lang2014"
numGenes <- filter(geneNumbers, Species == "Sac")$NumGenes
#note this is the "P1" dataset
names(Lang2014)[12:23] <- generations
#change the column names of P2 dataset
names(Lang2014)[24:35] <- paste0("P2_", generations)


#define empty vectors
numLineages <- c()
num_parallel_genes <- c()
num_non_parallel_genes <- c()
parallel_genes <- c()
c_hyper <- c()
p_chisq <- c()
estimate <- c()

#Loop over all available generations
for (g in generations){
  #create a new dataframe that filters the genes with a final proportion > 0 and only keep the columns we're intersted in for now
  print(g)
  Lang2014.g <- Lang2014 %>% 
    select(Population, Gene, Prop = g) %>% 
    filter(Prop >0) %>%
    filter(Gene != "Intergenic") %>% 
    arrange(Gene)
  
  #number of unique genes
  num_genes <- length((unique(Lang2014.g$Gene)))
  
  #number of lineages
  num_lineages <- length(unique(Lang2014.g$Population))
  numLineages <- append(numLineages, num_lineages)
  #set up matrix
  #add labels  onto the matrix ("Lineage" are columns, "Gene" are the row)
  data <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(Lang2014.g$Gene), unique(Lang2014.g$Population)))
  
  #subset the gene data (Lang2014) for each lineage in turn 
  #for each lineage, update the column adding 1's into the rows that correspond to genes that were identified in the dataset
  #which row corresponds to the genes in sub
  #turn each into a 1
  for(i in 1:num_lineages) {
    sub <- subset(Lang2014.g, Lang2014.g$Population == unique(Lang2014.g$Population)[i])
    sub2 <- subset(sub, Prop > 0)
    geneRows <- which(row.names(data) %in% sub2$Gene)
    data[geneRows, i] <- 1
    num_parallel <- data.frame(data, Count=rowSums(data, na.rm = FALSE, dims = 1), Gene = row.names(data))
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
  
  num_parallel_genes <- append(num_parallel_genes, nrow(genes_parallel))
  num_non_parallel_genes <- append(num_non_parallel_genes, nrow(Non_genes_parallel))
  parallel_genes <- append(parallel_genes, paste0(genes_parallel$Gene, collapse=", ")) 
  
  
  full_matrix <- rbind(data, array(0,c(numGenes-total_genes,ncol(data)))) 
  
  
  write_csv(Lang2014.g, path = paste0("data-out/Lang_", g, ".csv"))
  
  #run dgconstraint
  c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
  p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
  estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))
}

Lang2014_df <- tibble(paper ="Lang2014", environment= "YPD", as.numeric(generations), c_hyper, p_chisq, estimate, N.Populations=numLineages,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
write_csv(Lang2014_df, path = paste0("data-out/", paper, "_Analysis.csv"))

#plot the c_hyper by the number of generations 
#plot the number of parallel genes and the number of non-parallel genes on the same plot by the number of generations

Lang2014_graph1 <-Lang2014_df %>%
                    ggplot(aes(x=Generation,y=c_hyper)) + geom_point(size = 3) + theme(axis.text=element_text(size=13), axis.title=element_text(size=16,face="bold")) + scale_x_discrete(labels= X.label) + xlab("Time(generations)")


Parallel <-data.frame(Generation, Lang2014_df$N_genes.parallel) 
Non_Parallel <- data.frame(Generation, Lang2014_df$N_genes.notParallel)
Parallel.1 <- cbind(Parallel, rep("Parallel",12))
colnames(Parallel.1) <- c("Generations","Number_Genes","Gene_Type")
Non_Parallel.1 <- cbind(Non_Parallel, rep("Non_Parallel", 12))
colnames(Non_Parallel.1) <- c("Generations","Number_Genes","Gene_Type")
Data <- rbind(Parallel.1, Non_Parallel.1)

Lang2014_graph2 <- ggplot(Data, aes(x=Generations , y = Number_Genes , colour= Gene_Type)) +geom_point(size = 3) + theme(axis.text=element_text(size=13), axis.title=element_text(size=16,face="bold")) + scale_x_discrete(labels= X.label) + xlab("Time(generations)") + ylab("Number of Genes")

theme_set(theme_classic())

setwd("/Users/Margot/Desktop/Co-op/data-out/")
png("Lang2014_graph.1.png",height=10, width = 10, units="in", res=300)
Lang2014_graph.1
dev.off()
png("Lang2014_graph.2.png",height=10, width = 10, units="in", res=300)
Lang2014_graph.2
dev.off() 
