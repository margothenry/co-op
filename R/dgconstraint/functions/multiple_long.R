#' Calculations for Multiple Long Dataset
#' 
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a dataset with multiple selective pressures.
#' 
#' @param paper The name of the paper containing the dataset of interest.
#' @param dataset_name The actual name of the dataset (the part before "_usable.csv")
#' @param environment The environment in which the experiment occured.
#' @param generations Timepoint(s) in the data, if generations are used to notate. Must be numeric.
#' @param selective_pressure A list of the selective pressures in the data. i.e: temperatures, media, stressors.
#' @param species Specify if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species when prompted.
#' @param ploidy Haploid, diploid, etc. For E. coli, it's always haploid. 
#' (190620: In development) @param collapseMutations Specifies whether to run the analysis at the gene level or on distinct mutations within a gene. The default is at the gene level, i.e. to collapse all different mutations within a gene to one entry in the analysis.
#' @param numgenes The number of genes of the investigated species. If the species specified above is in the database, there's no need to enter a number here.
#' @param strain_info The specifics of the strain (i.e. the "mucoid" in "Wielgoss2016_mucoid")
#' @param days Timepoint(s) in the data, if days are used to notate. Remember to call with "days = ". Must be numeric.
#' @param flasks Timepoint(s) in the data, if flasks are are used to notate. Remember to call with "flasks = ". Must be numeric. Only 1 of the 3 potential timepoint types shall be called.
#' @param who_analyzed Who analyzed this dataset? Use initials of 1st & last names.
#' @return A table with all the calculated information.
#' @export 
#' @examples [update]
#'
multiple_long <- function(paper, dataset_name, environment, generations = NA, selective_pressure, species  = NA, ploidy, numgenes = NA, 
                          strain_info = NA, days = NA, flasks = NA, who_analyzed){

  geneNumbers <- read_csv(file.path(getwd(),"R/dgconstraint/inst/GeneDatabase.csv"), col_types = cols())

  data <- read_csv(file.path(getwd(), "data_in", "for_func", paste0(dataset_name, ".csv")), col_types = cols())

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
  c_hyper <- c()
  p_chisq <- c()
  estimate <- c()
  
  data.1 <- data %>% 
    arrange(gene) %>%
    drop_na(gene) %>%
    drop_na(population) %>%
    select(selective_pressure, population, gene, frequency)

  # The "for(j in selective_pressure)" loop is for multiple generations. The rest is essentially the same as single_long:
  cat("Evaluating constraint in ")
  for(j in selective_pressure) {
    cat("  ")
    cat(j)
    data.j <- data.1 %>% 
      filter(selective_pressure == j)
    
    num_genes <- length((unique(data.j$gene)))
    num_lineages <- length(unique(data.j$population))
    data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.j$gene), unique(data.j$population)))
  
  # In development.    
  # (TL): "collapseMutations" has a default value of TRUE:
  # if(collapseMutations){
  #   multiple_entry_genes <- subset(table(data.1$gene), table(data.1$gene) > 1)
  #   
  #   # These are our genes with only a single mutation:
  #   single_mutation_genes <- subset(data.1, gene %nin% names(multiple_entry_genes))  
  #   single_mutation_genes <- single_mutation_genes[, c(unique(population), "gene")]
  #   
  #   # These are the genes with multiple mutations. It may be the case in the future or in some circumstances that you want to know parallelism at the mutation rather than gene level. In that case don't include this.
  #   multiple_mutation_genes <- subset(data.1, gene %in% names(multiple_entry_genes))  
  #   multi_genes_matrix <- data.frame(gene = names(multiple_entry_genes))
  #   for (k in 1:length(multiple_entry_genes)){
  #     sub <- subset(data.1, gene == names(multiple_entry_genes)[k])
  #     for (j in unique(population)){
  #       # (TL): Fills "multi_genes_matrix" with k & j values from the loops. 
  #       # (TL): Each (k,j) co-ordinate is filled with the sum of all the values in the same population, thereby collapsing all potential cases of multiple mutations within the same gene.
  #       multi_genes_matrix[k, j] <- sum(sub[1:nrow(sub), j])
  #     }
  #   }
  #   # (TL) "multi_genes_matrix" is bound to "data.1". The multi-mutation genes are now treated as single-mutation genes.
  #   data.1 <- rbind(single_mutation_genes, multi_genes_matrix)
  #   # (TL): Rearrange "data.1" by gene names A-Z (to facilitate any changes from the newly-bound "multi_genes_matrix").
  #   data.1 <- data.1 %>% 
  #     arrange(gene) %>% 
  #     # (TL): Remove the columns that have nothing but NA values:
  #     filter(Reduce(`+`, lapply(., is.na)) != ncol(.))
  # }
  # 
  for(i in 1:num_lineages) {
    # (TL): Subset by population (via variable "sub").
    sub <- subset(data.j, data.j$population == unique(data.j$population)[i])
    # (TL): Eliminate those whose frequencies are 0 (via variable "sub2").
    sub2 <- subset(sub, frequency > 0)
    # (TL): In "data.array", assign to "geneRows" only the rows whose names are in "sub2".
    geneRows <- which(row.names(data.array) %in% sub2$gene)
    # (TL): These rows (for this particular iteration of "i") get the value of 1, essentially clumping all mutations into binary form (1 for presence):
    data.array[geneRows, i] <- 1
    # (TL): Creates a dataframe from the "population" column of "data.1", a "Count" column from all the values under said "population" column, and the row names are from the "gene" column in "data.1":
    num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), genes = row.names(data.array))
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
  # (TL): Get the values in the "genes" column, separated by ", ":
  parallel_genes <- append(parallel_genes, paste0(genes_parallel$genes, collapse=", ")) 
  # (TL): The full matrix consists of data.array and an all-0 array of (numgenes - total_genes) rows x (ncol(data.array)) columns:
  full_matrix <- rbind(data.array, array(0,c(numgenes-total_genes,ncol(data.array))))
  
  newdir <- file.path(getwd(), "data_out")
  if (!file.exists(newdir)){
    dir.create(newdir, showWarnings = FALSE)
    cat(paste("\n\tCreating new directory: ", newdir), sep="")
  }
  
  filename1 <- file.path(getwd(), "data_out", "intermediate", paste0("/", dataset_name, "_", j, ".csv"))
  write.csv(data.j, file=filename1, row.names=FALSE)
  
  c_hyper <- suppressWarnings(append(c_hyper, pairwise_c_hyper(full_matrix)))
  p_chisq <- suppressWarnings(append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200)))
  estimate <- suppressWarnings(append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = F)))
  
  c_hyper[c_hyper <= 0] <- 0
  c_hyper[c_hyper == "NaN"] <- 0
  }
  # (TL): The "func" column here will help with subsetting analyses for visualization:
  df <- tibble(paper = paper, dataset_name = dataset_name, environment = environment, generation = generations, day = days, flask = flasks, selective_pressure = selective_pressure, species = species, 
               ploidy = ploidy, strain_info = strain_info, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3), 
               N_genes.notParallel = num_non_parallel_genes, N_genes.parallel = num_parallel_genes, parallel_genes, who_analyzed, func = "multiple_long")
  filename2 <- file.path(getwd(), "data_out", "analyses", paste(dataset_name, "_Analysis.csv", sep=""))
  
  cat("\n")
  cat(paste("Writing file: ", filename2))
  write.csv(df, file=filename2, row.names=FALSE)
  cat("\n")
  }
  