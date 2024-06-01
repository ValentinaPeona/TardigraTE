### Small R script to create summary table of how many
### consensus sequences exist in a library per type of transposable element

### Usage: Rscript FigureLengthComparison.R <outfile> <list input_curated.lib> <list input_uncurated.lib>

### Valentina Peona, 12 Dec 2022

args = commandArgs(trailingOnly=TRUE)

# parse args
outfile = args[1]
files_curated = args[2:(((length(args)-1)/2)+1)]
files_uncurated = args[(((length(args)-1)/2)+2):length(args)]

library(ggplot2)
library(RColorBrewer)

table = data.frame(Consensus = character(), Curated_len = integer(), Uncurated_len = integer(), Species = character(), stringsAsFactors = FALSE)

for(i in 1:length(files_curated)){

	# the option of using # as comment char facilitates the comparison of names
	# keep it like that
	# the names after the # will not be imported into the table
	curated = read.table(files_curated[i], stringsAsFactors = FALSE, header = FALSE)
	uncurated = read.table(files_uncurated[i], stringsAsFactors = FALSE, header = FALSE)

	# remove the .raw label
	curated$V1 = sub(pattern = ".raw", replacement = "", x = curated$V1)
	uncurated$V1 = sub(pattern = ".raw", replacement = "", x = uncurated$V1)

	# further clean sequence names
	curated$V1 = sub(pattern = "_TTTT-TAT|_NA", replacement = "", x = curated$V1)
	curated$V1 = sub(pattern = "a$|b$", replacement = "", x = curated$V1)

	# re-format the curated and uncurated tables
	uncurated = data.frame(Name = uncurated$V1[seq(from = 1, to = nrow(uncurated), by = 2)], Sequence = uncurated$V1[seq(from = 2, to = nrow(uncurated), by = 2)], stringsAsFactors = FALSE)
	uncurated$Length = nchar(uncurated$Sequence)
	curated = data.frame(Name = curated$V1[seq(from = 1, to = nrow(curated), by = 2)], Sequence = curated$V1[seq(from = 2, to = nrow(curated), by = 2)], stringsAsFactors = FALSE)
	curated$Length = nchar(curated$Sequence)

	# table with lengths of the consensi for curates and uncurated libs


	if(grepl(pattern = "HD", x = files_curated[i])){

		uncurated$Name = sub(pattern = "nHd", replacement = "hypDuj", x = uncurated$Name)
		species = "HD"

	} else {

		species = "RV"

	}

	# summary
	summary = data.frame(Consensus = curated$Name, Curated_len = curated$Length, Uncurated_len = 0, Species = species, stringsAsFactors = FALSE)

	# get length for the uncurated version
	o = match(x = summary$Consensus, table = uncurated$Name)
	summary$Uncurated_len = uncurated$Length[o]

	table = rbind(table, summary)

}



	# get values for gradient for the scatterplot
	gradient = summary$Curated_len - summary$Uncurated_len

	# scatterplot
#scatterplot1 = ggplot(data = summary, aes(x = Uncurated_len, y = Curated_len)) + geom_point(aes(colour = gradient)) + scale_colour_gradient(low = "#7ccba2", high = "#045275") + theme_bw()
	scatterplot2 = ggplot(data = summary, aes(y = Uncurated_len, x = Curated_len)) + geom_point(aes(colour = gradient)) + scale_colour_gradient(low = "#7ccba2", high = "#045275") + theme_bw()

# remove the sequences that were not curated in both versions
boo = table$Curated_len == table$Uncurated_len
table = table[!boo,]
table = na.omit(table)

scatterplot1 = ggplot(data = table, aes(x = Uncurated_len, y = Curated_len)) + geom_point(aes(colour = Species)) + theme_bw() + xlim(NA, 11000) + ylim(NA, 11000)
ggsave(plot = scatterplot1, filename = outfile, units = "cm", device = "pdf", width = 20, height = 20, dpi = 300)

#f7feae,#b7e6a5,#7ccba2,#46aea0,#089099,#00718b,#045275