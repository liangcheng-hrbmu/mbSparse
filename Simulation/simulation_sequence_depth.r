library(HMP16SData)
V13_stool <-
  V13() %>%
  subset(select = HMP_BODY_SUBSITE == "Stool")
stool_16S <- V13_stool@assays@data[[1]]
dim(stool_16S)
sum(rowSums(stool_16S) != 0 )

taxa_map <- rowData(V13_stool)

taxa_map[1,7]
rownames(taxa_map)[1]

# Process the data and find the genus that overlaps with the Metagenomic data.

nz_prop <- apply(stool_16S, 1, FUN = function(x){
  sum(x == 0) / length(x)
})

set.seed(1234)
selected_taxa <- sample(which(nz_prop < 0.7), 300)

stool_16S <- stool_16S[selected_taxa,]

stool_16S <- t(stool_16S)

sum(colnames(stool_16S) %in% V13_stool@metadata$phylogeneticTree$tip.label)

phy_tree_taxa <- which(V13_stool@metadata$phylogeneticTree$tip.label %in% colnames(stool_16S))

stool_16S <- stool_16S[,V13_stool@metadata$phylogeneticTree$tip.label[phy_tree_taxa]]

selected_sample <- which(rowSums(stool_16S) >= 600 & rowSums(stool_16S) < 1500)

stool_16S <- stool_16S[selected_sample,]

total_reads <- rowSums(stool_16S)
total_reads <- (total_reads - min(total_reads))/sd(total_reads)

sequencing_depth <- 1000
stool_16S_normalized <- stool_16S / rowSums(stool_16S) * 10^3
stool_16S_processed <- log10(stool_16S_normalized+1)
# learn zero introduction parameter k from the benchmark truth.
zp = apply(stool_16S, 2, FUN = function(x){
  return( sum(x == 0)/length(x) )
})
mean_val = apply(stool_16S_processed, 2, FUN = function(x){
  return( mean(x) )
})
zp_map <- cbind(zp, mean_val)

# cor(log(zp), sqrt(nz_mean))

benchmark_truth <- stool_16S_normalized
for(j in 1:dim(stool_16S_normalized)[2]){
  abundance_truth <- log10(stool_16S_normalized[,j] + 1)
  # change the zero values to normal distribution
  abundance_truth[abundance_truth < log10(1.01)] = rnorm(sum(abundance_truth < log10(1.01)), mean = mean(abundance_truth[abundance_truth > log10(1.01)]), sd = sd(abundance_truth[abundance_truth > log10(1.01)]) ) 
  benchmark_truth[,j] = floor(10^abundance_truth - 1)
  benchmark_truth[benchmark_truth < 0] = 0
}

sum(benchmark_truth == 0) / (dim(benchmark_truth)[1] * dim(benchmark_truth)[2])
benchmark_truth_normalized <- log10( benchmark_truth / rowSums(benchmark_truth) * 10^3 + 1)

# adjust sequencing depth
run_simulation <- function(benchmark_truth, sequencing_depth, zp_map){
  simulated_truth <- benchmark_truth
  for(i in 1:dim(simulated_truth)[1]){
    simulated_truth[i,] = rmultinom(1, size = sequencing_depth, prob = benchmark_truth[i,] / sum(benchmark_truth[i,]))
  }
  simulated_truth <- floor(simulated_truth)
  #cat(rowSums(simulated_truth))
  simulated_truth_normalized <- simulated_truth / rowSums(simulated_truth) * 10^3
  simulated_truth_processed <- log10(simulated_truth_normalized+1)
  
  simulated_ZI_data <- simulated_truth_processed
  for(j in 1:dim(simulated_truth_processed)[2]){
    # introduce zeros to the simulated truth processed
    new_mean <- mean(simulated_truth_processed[,j])
    zp_idx <- which( zp_map[,2] + 0.5 > new_mean & zp_map[,2] - 0.5 < new_mean )
    zp_intro <- sample(zp_map[zp_idx,1], 1)
    new_abundance <-  rbinom(length(simulated_truth_processed[,j]), 1, prob = 1 - zp_intro) * simulated_ZI_data[,j]
    new_abundance[new_abundance == 0] = log10(1)
    simulated_ZI_data[,j] <- new_abundance
  }
  simulated_ZI_data_count = 10^simulated_ZI_data-1
  simulated_truth_processed_count = 10^simulated_truth_processed - 1
  full_normalize_str <- paste0("HMP16SData_", as.character(sequencing_depth), "_full_normalize.csv")
  zero_normalize_str <- paste0("HMP16SData_", as.character(sequencing_depth), "_zero_normalize.csv")
  zero_str <- paste0("HMP16SData_", as.character(sequencing_depth), "_zero.csv")
  write.csv(t(simulated_truth_processed), full_normalize_str)
  write.csv(t(simulated_ZI_data), zero_normalize_str)
  write.csv(t(simulated_ZI_data_count), zero_str)
}


set.seed(5)
run_simulation(benchmark_truth, 1000, zp_map)
run_simulation(benchmark_truth, 2000, zp_map)
run_simulation(benchmark_truth, 5000, zp_map)
run_simulation(benchmark_truth, 10000, zp_map)