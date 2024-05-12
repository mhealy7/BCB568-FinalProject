library(data.table)
library(car)
library(tidyverse)


#For the minimap AS probabilities
AS_data_pass <- fread("AS_results_pass.csv", header=TRUE, select = c("column", "average"))
AS_data_fail <- fread("AS_results_fail.csv", header=TRUE, select = c("column", "average"))

#setting the name of the columns to be more specific
setnames(AS_data_pass, "column", "TaxaID")
setnames(AS_data_fail, "column", "TaxaID")

setnames(AS_data_pass, "average", "average_pass")
setnames(AS_data_fail, "average", "average_fail")

#this is OK to do based on how the matrix is set up, so the number of rows (Taxa) match up and are ordered properly
AS_data <- merge(AS_data_pass, AS_data_fail, by="TaxaID")

As_data_no_taxa <- AS_data[,c("average_pass","average_fail")]

#have to pivot data to long
AS_data_long <- AS_data %>% 
  pivot_longer(
    cols = c(average_pass, average_fail), 
    names_to = "Condition",
    values_to = "Average"
)


leveneTest(Average ~ Condition, data = AS_data_long)
#variances are very equal


#Perform two-sample t-test
t.test(AS_data_pass$average_pass, AS_data_fail$average_fail, var.equal = TRUE)  #equal variances from levene test


#For Kraken

Kraken_report_pass <- fread("Kraken_report_pass.csv", header=TRUE, select = c("V6", "proportion"))
Kraken_report_fail <- fread("Kraken_report_fail.csv", header=TRUE, select = c("V6", "proportion"))

#setting the name of the columns to be more specific
setnames(Kraken_report_pass, "V6", "TaxaID")
setnames(Kraken_report_fail, "V6", "TaxaID")

setnames(Kraken_report_pass, "proportion", "proportion_pass")
setnames(Kraken_report_fail, "proportion", "proportion_fail")

#this is OK to do based since #rows is the same
kraken_data_merge <- merge(Kraken_report_pass, Kraken_report_fail, by="TaxaID")


#have to pivot data to long
kraken_data_long <- kraken_data_merge %>% 
  pivot_longer(
    cols = c(proportion_pass, proportion_fail), 
    names_to = "Condition",
    values_to = "Average"
  )

leveneTest(Average ~ Condition, data = kraken_data_long)
#variances are very equal


#Perform two-sample t-test
t.test(kraken_data_merge$proportion_pass, kraken_data_merge$proportion_fail, var.equal = TRUE)  #equal variances from levene test


