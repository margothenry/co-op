#' Calculations for Multiple Long Dataset
#' 
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a dataset with multiple selective pressures.
#' 
#' @param paper The data in .csv that you want to analyze.
#' @param environment The environment in which the experiment occured.
#' @param generation The generation the sequencing took place. Could also be a timepoint if the generation isn't specified in the paper. Make sure to include "units" (e.g. days, flasks) for non-generation entries.
#' @param selective_pressure A list of the selective pressures in the data. i.e: temperatures, media, stressors.
#' @param species Specify if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species when prompted.
#' @param numgenes The number of genes of the investigated species. If the species specified above is in the database, there's no need to enter a number here.
#' @return A table with all the calculated information.
#' @export 
#' @examples multiple_long("Jerison2017", c("YPD", "SC"), "500", c("OT", "HT"), "Sac")
#'
multiple_long <- function(paper, environment, generation, selective_pressure, species  = NA, numgenes = NA){

geneNumbers <- read_csv(file.path(getwd(),"dgconstraint/inst/GeneDatabase.csv"), col_types = cols())

data <- read_csv(file.path(getwd(), "data_in", "for_func", paste0(paper, ".csv")), col_types = cols())

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
  
  for(i in 1:num_lineages) {
    # (Tri): Subset by population (via variable "sub").
    sub <- subset(data.j, data.j$population == unique(data.j$population)[i])
    # (Tri): Eliminate those whose frequencies are 0 (via variable "sub2").
    sub2 <- subset(sub, frequency > 0)
    # (Tri): In "data.array", assign to "geneRows" only the rows whose names are in "sub2".
    geneRows <- which(row.names(data.array) %in% sub2$gene)
    # (Tri): These rows (for this particular iteration of "i") get the value of 1, essentially clumping all mutations into binary form (1 for presence):
    data.array[geneRows, i] <- 1
    # (Tri): Creates a dataframe from the "population" column of "data.1", a "Count" column from all the values under said "population" column, and the row names are from the "gene" column in "data.1":
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
  # (Tri): Get the values in the "genes" column, separated by ", ":
  parallel_genes <- append(parallel_genes, paste0(genes_parallel$genes, collapse=", ")) 
  # (Tri): The full matrix consists of data.array and an all-0 array of (numgenes - total_genes) rows x (ncol(data.array)) columns:
  full_matrix <- rbind(data.array, array(0,c(numgenes-total_genes,ncol(data.array))))
  
  newdir <- file.path(getwd(), "data_out")
  if (!file.exists(newdir)){
    dir.create(newdir, showWarnings = FALSE)
    cat(paste("\n\tCreating new directory: ", newdir), sep="")
  }
  
  filename1 <- file.path(getwd(), "data_out", "intermediate", paste0("/", paper, "_", j, ".csv"))
  write.csv(data.j, file=filename1, row.names=FALSE)
  
  c_hyper <- suppressWarnings(append(c_hyper, pairwise_c_hyper(full_matrix)))
  p_chisq <- suppressWarnings(append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200)))
  estimate <- suppressWarnings(append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = F)))
  
  c_hyper[c_hyper <= 0] <- 0
  c_hyper[c_hyper == "NaN"] <- 0
}
df <- tibble(paper = paper, environment = environment, generation = generation, selective_pressure, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3), 
             N_genes.notParallel = num_non_parallel_genes, N_genes.parallel = num_parallel_genes, parallel_genes)
filename2 <- file.path(getwd(), "data_out", "analyses", paste(paper, "_Analysis.csv", sep=""))
  
  cat("\n")
  cat(paste("Writing file: ", filename2))
  write.csv(df, file=filename2, row.names=FALSE)
  cat("\n")
  }
  