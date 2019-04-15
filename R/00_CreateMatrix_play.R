library(tidyverse)
#Read in the data
Lang2014 <- read_csv("data/Lang2014Nature.csv")

#create a new dataframe that filters the genes with a final proportion > 0 and only keep the columns we're intersted in for now
Lang2014.1000 <- Lang2014 %>% 
  filter(P1_1000 >0, Gene != "Intergenic") %>%
  select(Population, Gene, Prop = P1_1000) %>% 
  arrange(Gene)

#number of unique genes
num_genes <- length((unique(Lang2014.1000$Gene)))

#number of lineages
num_lineages <- length(unique(Lang2014.1000$Population))

#set up matrix
#add labels  onto the matrix ("Lineage" are columns, "Gene" are the row)
data <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(Lang2014.1000$Gene), unique(Lang2014.1000$Population)))
                                                                  
#subset the gene data (Lang2014) for each lineage in turn 
#for each lineage, update the column adding 1's into the rows that correspond to genes that were identified in the dataset
#which row corresponds to the genes in sub
#turn each into a 1
for(i in 1:num_lineages) {
  sub <- subset(Lang2014.1000, Lang2014.1000$Population == unique(Lang2014.1000$Population)[i])
  sub2 <- subset(sub, Prop > 0)
  geneRows <- which(row.names(data) %in% sub2$Gene)
  data[geneRows, i] <- 1
}

#write the matrix to the data folder ("Lang1000.csv")
write_csv(Lang2014.1000, path = "data/Lang1000.csv")

#Installing dgconstraint package
library(devtools)
install_github("samyeaman/dgconstraint")
library(dgconstraint)
pairwise_c_hyper(data)
allwise_p_chisq(data, num_permute = 100)
estimate_pa(data,ndigits = 4, show.plot = T)
