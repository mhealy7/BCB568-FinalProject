library(data.table)
library(Rsamtools)
library(ggplot2)
library(foreach)
library(doParallel)

path <- "/work/LAS/amytoth-lab/mhealy/Model_Building/06_MinimapAlign/SampleA/SampleA-barcode11-fail.bam"
bam <- scanBam(path)

bam_to_dt <- function(bamfile) {
  table <- data.table()
  
  param = ScanBamParam(tag = "AS",what=scanBamWhat())#telling R i actually want the data in my file
  bam <- scanBam(bamfile, param=param)#why is this program designed like this?
  
  
  table[, sname := bam[[1]]$qname]
  #table[, flag := bam[[1]]$flag]
  rname <- gsub(".*\\|(.*?)\\|.*", "\\1", bam[[1]]$rname) #extracts just taxaID
  table[, taxa := rname]
  #table[, strand := bam[[1]]$strand]
  #table[, pos := bam[[1]]$pos]
  #table[, qwidth := bam[[1]]$qwidth]
  table[, mapq := bam[[1]]$mapq]
  #table[, cigar := bam[[1]]$cigar]
  #table[, mrnm := bam[[1]]$mrnm]
  table[, AS := bam[[1]]$tag$AS]
  table <- unique(table)#really should figure out why my sam files are duplicated!
  return(table)
}


dt <- bam_to_dt(path)

#removing no matches
dt <- dt[!is.na(taxa),]

#keeping only the top 5 AS scores per read
new_dt <- dt[order(sname, -AS), head(.SD, 5), by = sname]

#check to verify
read_counts <- new_dt[, .(Count = .N), by = sname]

#replacing the old dt once its been confirmed new_dt is legit
dt <- new_dt

#creating the matrix based off the data
unique_sequence <- unique(dt$sname)
unique_taxa <- unique(dt$taxa)
my_matrix <- matrix(0, nrow = length(unique_sequence), ncol = length(unique_taxa),
                    dimnames = list(unique_sequence, unique_taxa))

#make sure to change 16 in case you do not (or maybe have more) have enough cores
numCores <- 16 - 1  #can only allocate 16 cores for an interactive Rstudio job on nova. -1 just to not overload anything
cl <- makeCluster(numCores)
registerDoParallel(cl)

results <- foreach(i = 1:nrow(my_matrix), .combine='rbind', .packages='data.table') %dopar% {
  sname_index <- rownames(my_matrix)[i]
  
  #Get the entries for the given read
  read_data <- dt[sname == sname_index]
  
  #get the max AS score for a given read
  max_AS <- max(read_data$AS, na.rm = TRUE)
  read_data[, AS := exp(AS - max_AS)]  #Apply softmax transformation replacing the old AS values for optimization
  
  # Create a vector to store results for this row
  row_result <- numeric(ncol(my_matrix))
  
  #looping through the taxa
  for (j in 1:ncol(my_matrix)) {
    taxa_index <- colnames(my_matrix)[j]
    max_score <- read_data[taxa == taxa_index, max(AS, na.rm = TRUE)]
    sum_score <- sum(read_data$AS)  #sum_score for the softmax'd AS
    
    #just incase any NAs
    max_score <- ifelse(is.na(max_score), 0, max_score)
    row_result[j] <- max_score / sum_score  
  }
  
  return(row_result)
}

# Stop the parallel cluster
stopCluster(cl)

# Assign results back to my_matrix
for (i in 1:nrow(my_matrix)) {
  my_matrix[i, ] <- results[i, ]
}

#Probability 0s (where a taxa does not occur for a read) result in -inf, i could add a check earlier in this stage but this is easy!
my_matrix[is.infinite(my_matrix) & my_matrix < 0] <- 0

column_averages <- colMeans(my_matrix)

#putting it into a table
results_table <- data.frame(column = colnames(my_matrix), average = column_averages)

#save as a RDS, to make it easily loadable to R later
saveRDS(my_matrix, "AS_matrix_fail.rds")

write.csv(results_table, "AS_results_fail.csv", row.names = TRUE)







