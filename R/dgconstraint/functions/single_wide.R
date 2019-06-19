#' Calculations for a Single Wide Dataset
#'
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a dataset with a single generation.
#' 
#' @param paper The data in .csv that you want to analyze.
#' @param environment The environment in which the experiment occured.
#' @param generations Timepoint(s) in the data, if generations are used to notate.
#' @param selective_pressure A list of the selective pressures in the data. i.e: temperatures, media, stressors.
#' @param species Specifies if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species when prompted.
#' @param ploidy Haploid, diploid, etc. For E. coli, it's always haploid. 
#' @param collapseMutations Specifies whether to run the analysis at the gene level or on distinct mutations within a gene. The default is at the gene level, i.e. to collapse all different mutations within a gene to one entry in the analysis.
#' @param numgenes The number of genes of the investigated species. If the species specified above is in the database, there's no need to enter a number here.
#' @param days Timepoint(s) in the data, if days are used to notate. Remember to call with "days = ".
#' @param flasks Timepoint(s) in the data, if flasks are are used to notate. Remember to call with "flasks = ". Only 1 of the 3 potential timepoint types shall be called.
#' @return A table with all the calculated information.
#' @export 
#' @examples: [update]
#'
single_wide <- function(paper, environment, generation = NA, selective_pressure, species = NA, ploidy, population, collapseMutations = TRUE, numgenes = NA, days = NA, flasks = NA){

  geneNumbers <- read_csv(file.path(getwd(),"R/dgconstraint/inst/GeneDatabase.csv"), col_types = cols())
  
  data <- read_csv(file.path(getwd(), "data_in", "for_func", paste0(paper, ".csv")), col_types = cols())
  
# (Tri): If the species name is found in the database, the "numgenes" can be retrieved there.
  if (species %in% geneNumbers$Species){
    numgenes <- filter(geneNumbers, Species == species)$NumGenes  
  }
  
# (Tri): Otherwise, enter "numgenes" by hand.
  if(is.na(numgenes)){
    prompt <- "Your species is unspecified or not in our database. How many genes does it have? \n"
    numgenes <- as.numeric(readline(prompt))
  }
  
#(Tri): Arrange gene names A-Z, and omit rows containing missing values in "gene" and/or "population":
  data.1 <- data %>%
    arrange(gene) %>%
    drop_na(gene) %>%
    drop_na(population)
  
  num_lineages <- length(unique(population))
  num_genes <- length((unique(data.1$gene)))
  
  data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.1$gene), unique(population)))
  
  # (Tri): "collapseMutations" has a default value of TRUE:
  if(collapseMutations){
    multiple_entry_genes <- subset(table(data.1$gene), table(data.1$gene) > 1)
  
    # These are our genes with only a single mutation:
    single_mutation_genes <- subset(data.1, gene %nin% names(multiple_entry_genes))  
    single_mutation_genes <- single_mutation_genes[, c(unique(population), "gene")]
    
    # These are the genes with multiple mutations. It may be the case in the future or in some circumstances that you want to know parallelism at the mutation rather than gene level. In that case don't include this.
    multiple_mutation_genes <- subset(data.1, gene %in% names(multiple_entry_genes))  
    multi_genes_matrix <- data.frame(gene = names(multiple_entry_genes))
    for (k in 1:length(multiple_entry_genes)){
      sub <- subset(data.1, gene == names(multiple_entry_genes)[k])
      for (j in unique(population)){
    # (Tri): Fills "multi_genes_matrix" with k & j values from the loops. 
    # (Tri): Each (k,j) co-ordinate is filled with the sum of all the values in the same population, thereby collapsing all potential cases of multiple mutations within the same gene.
        multi_genes_matrix[k, j] <- sum(sub[1:nrow(sub), j])
      }
    }
    # (Tri) "multi_genes_matrix" is bound to "data.1". The multi-mutation genes are now treated as single-mutation genes.
    data.1 <- rbind(single_mutation_genes, multi_genes_matrix)
    # (Tri): Rearrange "data.1" by gene names A-Z (to facilitate any changes from the newly-bound "multi_genes_matrix").
    data.1 <- data.1 %>% 
      arrange(gene) %>% 
    # (Tri): Remove the columns that have nothing but NA values:
       filter(Reduce(`+`, lapply(., is.na)) != ncol(.))
  }
  
  #this is from single_long
  #num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), genes = row.names(data.array))
  
  # (Tri): Creates a dataframe from the "population" column of "data.1", a "Count" column from all the values under said "population" column, and the row names are from the "gene" column in "data.1":
  # (Tri): Potential problem: The "gene" column might contain duplicates, which can be troublesome if we want to make row names out of them.
  num_parallel <- data.frame(data.1[, population], Count=rowSums(data.1[, population], na.rm = TRUE, dims = 1), row.names = data.1$gene, genes = row.names(data.array))
  
  genes_parallel <- num_parallel %>%
    as_tibble() %>%
    filter(Count > 1)
  
  non_genes_parallel <- num_parallel %>%
    as_tibble() %>%
    filter(Count == 1)
  
  num_parallel_genes <- nrow(genes_parallel)
  num_non_parallel_genes <- nrow(non_genes_parallel)
  total_genes <- num_non_parallel_genes + num_parallel_genes
  parallel_genes <- paste0(genes_parallel$genes, collapse=", ")
  
  extragenes <- data.frame(array(0, dim = c(numgenes - total_genes, length(population))))
  names(extragenes) <- names(num_parallel[, population])
  full_matrix <- rbind(num_parallel[, population], extragenes)
  
  c_hyper <- c()
  p_chisq <- c()
  estimate <- c()
  c_hyper <- suppressWarnings(append(c_hyper, pairwise_c_hyper(full_matrix)))
  p_chisq <- suppressWarnings(append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200)))
  estimate <- suppressWarnings(append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = T)))
  
  c_hyper[c_hyper <= 0] <- 0
  c_hyper[c_hyper == "NaN"] <- 0
  
  df <- tibble(paper = paper, environment = environment, generation = generation, day = days, flask = flasks, selective_pressure = selective_pressure, species = species, ploidy, 
               c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3), 
               N_genes.notParallel = num_non_parallel_genes, N_genes.parallel = num_parallel_genes, parallel_genes)
  
  newdir <- file.path(getwd(), "data_out", "intermediate")
  if (!file.exists(newdir)){
    dir.create(newdir, showWarnings = FALSE)
    cat(paste("\n\tCreating new directory: ", newdir), sep="")
  }
  
  filename <- file.path(getwd(), "data_out", "analyses", paste(paper, "_Analysis.csv", sep=""))
  cat("\n")
  cat(paste("Writing file: ", filename))
  write.csv(df, file = filename, row.names=FALSE)
  cat("\n")
}
