library(data.table)
library(ggplot2)


path <- "/work/LAS/amytoth-lab/mhealy/Model_Building/05_TaxaByBarcode/SampleA/"

file <- "SampleA-barcode11-pass-report.txt"

complete_path <- paste0(path,file)

report <- fread(complete_path)

V1 <- "Percentage of fragments covered by the clade rooted at this taxon"
V2 <- "Number of fragments covered by the clade rooted at this taxon"
V3 <- "Number of fragments assigned directly to this taxon"
V4 <- "A rank code"
V5 <- "NCBI taxonomic ID number"
V6 <- "Indented scientific name"

sum(report$V3)

#removes non terminal taxa
filtered_report <- report[V3 != 0]

#removes unmatched reads
filtered_report <- filtered_report[V6 != "unclassified"]

#removing non species reads
filtered_report <- filtered_report[V4 == "S"]

#get k taxa
k <- length(unique(filtered_report$V6))

k

#get n
n <- sum(filtered_report$V3)

n

#get proportions
filtered_report[, proportion := (V3 / n)]



#sum of probabilities
sum(filtered_report$proportion)

saveRDS(filtered_report, "Kraken_report_pass.rds")
write.csv(filtered_report, "Kraken_report_pass.csv")


