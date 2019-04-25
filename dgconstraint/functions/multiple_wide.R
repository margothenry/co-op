#' Calculations for a Multiple Wide Dataset 
#'
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a data set with multiple generations.
#' 
#' @param paper the data in csv that you want to analyze, in a folder named data_in
#' @param environment The environment in which the experiment occured
#' @param species Specify if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species
#' @param generations a list of generations in the data
#' @return a table with all the calculated infromation
#' @export 
#' @examples 
#' multiple_wide("Author2018","YPD", "Sac", c("0", "100", "500" , "1000"))
#' 
multiple_wide <- function(paper, generations, environment, species = NA, numgenes = NA){

geneNumbers <- read_csv(file.path(getwd(),"dgconstraint/inst/geneDatabase.csv"), col_types = cols())

data <- read_csv(file.path(getwd(), "data_in", paste0(paper, ".csv")), col_types = cols())

if (species %in% geneNumbers$Species){
  numgenes <- filter(geneNumbers, Species == species)$NumGenes  
}

if(is.na(numgenes)){
  prompt <- "Your species is unspecified or not in our database. How many genes does it have? \n"
  numgenes <- as.numeric(readline(prompt))
}
  numLineages <- c()
  num_parallel_genes <- c()
  num_non_parallel_genes <- c()
  parallel_genes <- c()
  c_hyper <-
    c()
  p_chisq <- c()
  estimate <- c()
  
  cat("Evaluating constraint in ")
  for (g in generations) {
    cat("  ")
    cat(g)
    data.1 <- data %>% 
      arrange(gene) %>%
      drop_na(gene) %>%
      drop_na(population)
    
    data.g <- data.1 %>% 
      select(population, gene, prop = g) %>% 
      filter(prop >0)
    
    #number of unique genes
    num_genes <- length((unique(data.g$gene)))
    
    #number of lineages
    num_lineages <- length(unique(data.g$population))
    numLineages <- append(numLineages, num_lineages)
    
    #making the matrix
    data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.g$gene), unique(data.g$population)))
    
    for(i in 1:num_lineages) {
      sub <- subset(data.g, data.g$population == unique(data.g$population)[i])
      sub2 <- subset(sub, prop > 0)
      geneRows <- which(row.names(data.array) %in% sub2$gene)
      data.array[geneRows, i] <- 1
      num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), genes = row.names(data.array))
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
    parallel_genes <- append(parallel_genes, paste0(genes_parallel$genes, collapse=", ")) 
    
    
    full_matrix <- rbind(data.array, array(0,c(numgenes-total_genes,ncol(data.array))))
   
    newdir <- file.path(getwd(), "data_out")
    if (!file.exists(newdir)){
      dir.create(newdir, showWarnings = FALSE)
      cat(paste("\n\tCreating new directory: ", newdir), sep="")
    }
    
    filename1 <- file.path(getwd(), "data_out", paste0("/", paper, "_", g, ".csv"))
    write.csv(data.g, file=filename1, row.names=FALSE)
    
    c_hyper <- suppressWarnings(append(c_hyper, pairwise_c_hyper(full_matrix)))
    p_chisq <- suppressWarnings(append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200)))
    estimate <- suppressWarnings(append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T)))
    
    c_hyper[c_hyper <= 0] <- 0
    c_hyper[c_hyper == "NaN"] <- 0
  }
  
  df <- tibble( paper = paper, environment = environment, generations, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3) ,N_genes.notParallel= num_non_parallel_genes, N_genes.parallel=num_parallel_genes, parallel_genes)
  
  filename2 <- file.path(getwd(), "data_out", paste(paper, "_Analysis.csv", sep=""))
  cat("\n")
  cat(paste("Writing file: ", filename2))
  write.csv(df, file=filename2, row.names=FALSE)
  cat("\n")
}
