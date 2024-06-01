### Small R script to create summary table of how many
### consensus sequences exist in a library per type of transposable element

### Usage: Rscript TableOverviewConsensi.R outfile TE.lib

### Valentina Peona, 12 Dec 2022

args = commandArgs(trailingOnly=TRUE)

outfile = args[1]
files = args[2:length(args)]

library(dplyr)

# read library

for(filename in files){

 library = read.table(file = filename, stringsAsFactors = FALSE, header = FALSE, comment.char = "")

 # get only fasta headers
 boo = grepl(pattern = ">", x = library$V1)
 library = library$V1[boo]

 # remove '>' character
 library = sub(pattern = ">", x = library, replacement = "")

 # remove what is not a transposable element
 boo = grepl(pattern = "ARTEFACT|Low_complexity|Simple_repeat|Satellite|tRNA|rRNA", x = library)
 library = library[!boo]

 # get type of TE (classification before the "/")
 types = data.frame(Values = sapply(strsplit(x = sapply(strsplit(x = library, split = "/"), "[", 1), split = "\\#"), "[", 2))

 # remove the '?'
 types$Values = sub(pattern = "?", x = types$Values, replacement = "", fixed = TRUE)

 # create summary
 summary = types %>% group_by(Values) %>% tally()

 # add first row with info about the library
 sink(file = outfile, append = TRUE)
 cat(paste("\n", filename, "\n"))
 sink()

 # write results
 write.table(x = summary, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

}