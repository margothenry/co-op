#' Calculations for a Single Wide Dataset
#'
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a data set with a single generation.
#' 
#' @param paper the data in csv that you want to analyze, in a folder named data_in
#' @param environment The environment in which the experiment occured
#' @param species Specify if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species
#' @param population a list of populations in the dataset
#' #' @param collapseMutations specifys whether to run analysis at the level of the gene or on distinct mutations within a gene. The default is at the gene level, i.e. to collapse all different mutations within a gene to one entry in the analysis.
#' @return a table with all the calculated infromation
#' @export 
#' @examples 
#' single_wide("Author2018","YPD", "Sac", c("P1", "P2", "P3" ,"P4", "P5"))

single_wide <- function(paper, environment, population, species = NA, collapseMutations = TRUE, numgenes = NA){
  geneNumbers <- read_csv(file.path(getwd(),"inst/geneDatabase.csv"), col_types = cols())
  
  data <- read_csv(file.path(getwd(), "data_in/original", paste0(paper, ".csv")), col_types = cols())
  
  if (species %in% geneNumbers$Species){
    numgenes <- filter(geneNumbers, Species == species)$Numgenes  
  }
  
  if(is.na(numgenes)){
    prompt <- "Your species is unspecified or not in our database. How many genes does it have? \n"
    numgenes <- as.numeric(readline(prompt))
  }
  
  data.1 <- data %>%
    arrange(gene) %>%
    drop_na(gene) %>%
    drop_na(population)
  
  num_lineages <- length(unique(population))
  num_genes <- length((unique(data.1$gene)))
  
  #I don't think this is used atm
  data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.1$gene), unique(population)))
  
  if(collapseMutations){
    multiple_entry_genes <- subset(table(data.1$gene), table(data.1$gene) >1)
    
    # These are our genes with ony a single mutation
    single_mutation_genes <- subset(data.1, gene %nin% names(multiple_entry_genes))  
    single_mutation_genes <- single_mutation_genes[, c(unique(population), "gene")]
    
    # These are the genes with multiple mutations. It may be the case in the future or in some circumstances that you want to know parallelism at the mutation rather than gene level. In that case don't include this.
    multiple_mutation_genes <- subset(data.1, gene %in% names(multiple_entry_genes))  
    multi_genes_matrix <- data.frame(gene = names(multiple_entry_genes))
    for (k in 1:length(multiple_entry_genes)){
      sub <- subset(data.1, gene == names(multiple_entry_genes)[k])
      for (j in unique(population)){
        multi_genes_matrix[k, j] <- sum(sub[1:nrow(sub), j])
      }
    }
    
    data.1 <- rbind(single_mutation_genes, multi_genes_matrix)
    data.1 <- data.1 %>% 
      arrange(gene) 
  }
  
  
  num_parallel <- data.frame(data.1[, population], Count=rowSums(data.1[, population], na.rm = FALSE, dims = 1),row.names= data.1$gene)
  
  genes_parallel <- num_parallel %>%
    as_tibble() %>%
    filter(Count > 1)
  
  Non_genes_parallel <- num_parallel %>%
    as_tibble() %>%
    filter(Count == 1)
  
  num_parallel_genes <- nrow(genes_parallel)
  num_non_parallel_genes <- nrow(Non_genes_parallel)
  total_genes <- num_non_parallel_genes + num_parallel_genes
  
  extragenes <- data.frame(array(0, dim = c(numgenes-total_genes, length(population))))
  names(extragenes) <- names(num_parallel[, population])
  full_matrix <- rbind(num_parallel[, population], extragenes)
  
  c_hyper <- c()
  p_chisq <- c()
  estimate <- c()
  c_hyper <- append(c_hyper, pairwise_c_hyper(full_matrix))
  p_chisq <- append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200))
  estimate <- append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T))
  
c_hyper[c_hyper <= 0] <- 0
c_hyper[c_hyper == "NaN"] <- 0

df <- tibble( paper = paper, environment = environment, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)

newdir <- file.path(getwd(), "data_out")
if (!file.exists(newdir)){
  dir.create(newdir, showWarnings = FALSE)
  cat(paste("\n\tCreating new directory: ", newdir), sep="")
}

filename <- file.path(getwd(), "data_out", paste(paper, "_Analysis.csv", sep=""))
write.csv(df, file=filename, row.names=FALSE)

}
