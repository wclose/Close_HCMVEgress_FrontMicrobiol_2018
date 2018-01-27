library(tidyverse)

# import the sam tables for each comparison
# tvm = towne vs mock, avm = ad169 vs mock, 12/96 = hpi
# also import the data file from the microarray
tvm12 <- read.table("tvm12", header = T, sep = "\t")
tvm96 <- read.table("tvm96", header = T, sep = "\t")
avm12 <- read.table("avm12", header = T, sep = "\t")
avm96 <- read.table("avm96", header = T, sep = "\t")
pellett_probe <- read.table("pellett_probe", header = T, sep = "\t")

# change data from factors to chrs for join
tvm12$PROBE_ID <- as.character(tvm12$PROBE_ID)
tvm96$PROBE_ID <- as.character(tvm96$PROBE_ID)
avm12$PROBE_ID <- as.character(avm12$PROBE_ID)
avm96$PROBE_ID <- as.character(avm96$PROBE_ID)
pellett_probe$PROBE_ID <- as.character(pellett_probe$PROBE_ID)

tvm12$Symbol <- as.character(tvm12$Symbol)
tvm96$Symbol <- as.character(tvm96$Symbol)
avm12$Symbol <- as.character(avm12$Symbol)
avm96$Symbol <- as.character(avm96$Symbol)
pellett_probe$SYMBOL <- as.character(pellett_probe$SYMBOL)

# create lists of shared probe names
tvm_sam_probe_list <- full_join(tvm12, tvm96)
avm_sam_probe_list <- full_join(avm12, avm96)
sam_probe_list <- full_join(avm_sam_probe_list, tvm_sam_probe_list)

# join the lists of probe names to the data for each
tvm_sam_probe_data <- left_join(tvm_sam_probe_list, pellett_probe, by = "PROBE_ID")
avm_sam_probe_data <- left_join(avm_sam_probe_list, pellett_probe, by = "PROBE_ID")
sam_probe_data <- left_join(sam_probe_list, pellett_probe, by = "PROBE_ID")

# change symbol col name from "Symbol" to "alt_Symbol" because there are slight differences
# in names used and the "SYMBOL" col has better annotation
colnames(sam_probe_data)[colnames(sam_probe_data) == "Symbol"] <- "alt_Symbol"

# save the resulting lists and datasets as .txt files for later use
write.table(tvm_sam_probe_list, "tvm_sam_probe_list.txt", sep = "\t")
write.table(tvm_sam_probe_data, "tvm_sam_probe_data.txt", sep = "\t")

write.table(avm_sam_probe_list, "avm_sam_probe_list.txt", sep = "\t")
write.table(avm_sam_probe_data, "avm_sam_probe_data.txt", sep = "\t")

write.table(sam_probe_list, "sam_probe_list.txt", sep = "\t")
write.table(sam_probe_data, "sam_probe_data.txt", sep = "\t")

# ready to calculate the dellog2fc for the data (uses del_log2_fc.R)