### Small R script to create summary table of how many
### consensus sequences exist in a library per type of transposable element

### Usage: Rscript TablePutativeActive.R outfile.tsv rmsk.out 

### Valentina Peona, 20 January 2023

# ten intact copies at 0-1%

args = commandArgs(trailingOnly=TRUE)

outfile = args[1]
filenames = args[2:length(args)]


library(data.table)
library(dplyr)

for(file in filenames){

print(file)

data = fread(file = file, header = F, stringsAsFactors = F, skip = 3, fill = T)
data = data[,c(5,6,7,8,9,10,11,2,15)]
names(data) = c("Scaffold", "Begin", "End", "Remain", "Strand", "Element", "Family", "Divergence", "ID")

# group elements in major classes
data$newFamily = sapply(strsplit(x = data$Family,split = "/"), "[", 1)

# keep only these elements
KEEP = c("SINE", "LINE", "LTR", "DNA", "Unknown", "RC")
data = data[data$newFamily %in% KEEP,] 

# remove parentheses from the Remain column
data$Begin = as.integer(gsub(pattern = "[()]", x = as.character(data$Begin), replacement = ""))
data$End= as.integer(gsub(pattern = "[()]", x = as.character(data$End), replacement = ""))
data$Remain = as.integer(gsub(pattern = "[()]", x = as.character(data$Remain), replacement = ""))

# calculate the length of the elements
data$tot_len = abs(data$Begin - data$End) + data$Remain

# calculate coverage of the element
data$coverage = (abs(data$Begin - data$End) / data$tot_len) * 100

# get consensi with copies at less than 1%
data_filtered = data[data$Divergence < 1,]
summary = data_filtered %>% group_by(Element) %>% tally() %>% arrange(desc(n)) %>% filter(n >= 10)

sink(file = outfile, append = TRUE)
print(cat("\n", file, "\n"))
sink()
write.table(x = summary, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

}