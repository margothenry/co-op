#' Calculations for a Single Long Dataset
#'
#' This function allows you to calculate the pairwise C-score using the hypergeometric approach, a p-value for 'all lineages' contrast using chi-square, and the estimates of the effective proportion of adaptive loci for a data set with a single generation.
#' 
#' @param paper the data in csv that you want to analyze, in a folder named data_in
#' @param environment The environment in which the experiment occured
#' @param species Specify if the organism is "Sac" or "Ecoli_K12" or "Ecoli_O157-H7", or manually input the gene count of your species
#' @return a table with all the calculated infromation
#' @export 
#' @examples 
#'single_long("Author2018","YPD", "Sac")
#'
single_long <- function(paper, environment, generation, selective_pressure, species = NA, numgenes = NA){
  
  geneNumbers <- read_csv(file.path(getwd(),"dgconstraint/inst/GeneDatabase.csv"), col_types = cols())
  data <- read_csv(file.path(getwd(), "data_in", "for_func", paste0(paper, ".csv")), col_types = cols())   

if (species %in% geneNumbers$Species){
  numgenes <- filter(geneNumbers, Species == species)$NumGenes  
}

if(is.na(numgenes)){
  prompt <- "Your species is unspecified or not in our database. How many genes does it have? \n"
  numgenes <- as.numeric(readline(prompt))
}

  data.1 <- data %>%
  arrange(gene) %>%
  drop_na(gene) %>%
  drop_na(population)

num_genes <- length((unique(data.1$gene)))
num_lineages <- length(unique(data.1$population))

# (Tri): Create an array for genes in different lineages: 
# (Tri): This essentially converts the data from long to wide by making column names out of population names:
data.array <- array(0, dim =c(num_genes, num_lineages), dimnames = list(unique(data.1$gene), unique(data.1$population)))

for(i in 1:num_lineages) {
  # (Tri): Subset by population (via "sub").
  sub <- subset(data.1, data.1$population == unique(data.1$population)[i])
  # (Tri): Eliminate those whose frequencies are 0 (via "sub2").
  sub2 <- subset(sub, frequency > 0)
  # (Tri): In "data.array", assign to "geneRows" only the rows whose names are in "sub2".
  geneRows <- which(row.names(data.array) %in% sub2$gene)
  # (Tri): These rows (for this particular iteration of "i") get the value of 1:
  data.array[geneRows, i] <- 1
  # (Tri): Creates a dataframe from the "population" column of "data.1", a "Count" column from all the values under said "population" column, and the row names are from the "gene" column in "data.1":
  num_parallel <- data.frame(data.array, Count=rowSums(data.array, na.rm = FALSE, dims = 1), genes = row.names(data.array))
}

genes_parallel <- num_parallel %>%
  as_tibble() %>%
  filter(Count > 1)


Non_genes_parallel <- num_parallel %>%
  as_tibble() %>%
  filter(Count == 1)

num_parallel_genes <- nrow(genes_parallel)
num_non_parallel_genes <- nrow(Non_genes_parallel)
total_genes <- num_non_parallel_genes + num_parallel_genes
# (Tri): Get the values in the "genes" column, separated by ", ":
parallel_genes <- paste0(genes_parallel$genes, collapse=", ")

full_matrix <- rbind(data.array, array(0,c(numgenes-total_genes,ncol(data.array))))

c_hyper <- c()
p_chisq <- c()
estimate <- c()
c_hyper <- suppressWarnings(append(c_hyper, pairwise_c_hyper(full_matrix)))
p_chisq <- suppressWarnings(append(p_chisq, allwise_p_chisq(full_matrix, num_permute = 200)))
estimate <- suppressWarnings(append(estimate, estimate_pa(full_matrix,ndigits = 4, show.plot = F)))

c_hyper[c_hyper <= 0] <- 0
c_hyper[c_hyper == "NaN"] <- 0


df <- tibble(paper = paper, environment = environment, generation = generation, selective_pressure = selective_pressure, c_hyper = round(c_hyper, 3), p_chisq, estimate = round(estimate, 3), 
             N_genes.notParallel = num_non_parallel_genes, N_genes.parallel = num_parallel_genes, parallel_genes)

newdir <- file.path(getwd(), "data_out", "intermediate")
if (!file.exists(newdir)){
  dir.create(newdir, showWarnings = FALSE)
  cat(paste("\n\tCreating new directory: ", newdir), sep="")
}

filename <- file.path(getwd(), "data_out", "analyses", paste(paper, "_Analysis.csv", sep=""))
cat("\n")
cat(paste("Writing file: ", filename))
write.csv(df, file=filename, row.names=FALSE)
cat("\n")
}
