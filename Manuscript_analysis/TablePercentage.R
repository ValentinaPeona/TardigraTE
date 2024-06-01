### Small R script to create summary table of how many
### consensus sequences exist in a library per type of transposable element

### Usage: Rscript TablePercentage.R outfile TE.lib

### Valentina Peona, 20 January 2023

args = commandArgs(trailingOnly=TRUE)

outfile = args[1]
files = args[2:length(args)]

table = data.frame(Species = character(), Library = character(), LINE = character(), LINE_perc = character(), SINE = character(), SINE_perc = character(), LTR = character(), LTR_perc = character(), DNA = character(), DNA_perc = character(), Unknown = character(), Unknown_perc = character(), Total = character(), Total_perc = character())

for (i in 1:length(files)){

	species = strsplit(x = sub(pattern = ".tbl", replacement = "", x = files[i]), split = "_")[[1]][1]
	species = sub(pattern = "Intermediate/RMSK/", replacement = "", x = species)
	library = strsplit(x = sub(pattern = ".tbl", replacement = "", x = files[i]), split = "_")[[1]][2]

	# get absolute bps
	lines = system(command = paste0("grep 'LINEs:' ", files[i], " | tr -s ' ' | cut -f4 -d ' '"), intern = TRUE)
	sines = system(command = paste0("grep 'SINEs:' ", files[i], " | tr -s ' ' | cut -f4 -d ' '"), intern = TRUE)
	ltr = system(command = paste0("grep 'LTR elements:' ", files[i], " | tr -s ' ' | cut -f5 -d ' '"), intern = TRUE)
	dna = system(command = paste0("grep 'DNA transposons' ", files[i], " | tr -s ' ' | cut -f4 -d ' '"), intern = TRUE)
	unknown = system(command = paste0("grep 'Unclassified:' ", files[i], " | tr -s ' ' | cut -f3 -d ' '"), intern = TRUE)
	total = system(command = paste0("grep 'Total interspersed repeats:' ", files[i], " | tr -s ' ' | cut -f4 -d ' '"), intern = TRUE)

	# get percentages
	lines_perc = system(command = paste0("grep 'LINEs:' ", files[i], " | tr -s ' ' | cut -f6 -d ' '"), intern = TRUE)
	sines_perc = system(command = paste0("grep 'SINEs:' ", files[i], " | tr -s ' ' | cut -f6 -d ' '"), intern = TRUE)
	ltr_perc = system(command = paste0("grep 'LTR elements:' ", files[i], " | tr -s ' ' | cut -f7 -d ' '"), intern = TRUE)
	dna_perc = system(command = paste0("grep 'DNA transposons' ", files[i], " | tr -s ' ' | cut -f6 -d ' '"), intern = TRUE)
	unknown_perc = system(command = paste0("grep 'Unclassified:' ", files[i], " | tr -s ' ' | cut -f5 -d ' '"), intern = TRUE)
	total_perc = system(command = paste0("grep 'Total interspersed repeats:' ", files[i], " | tr -s ' ' | cut -f6 -d ' '"), intern = TRUE)

	# create row for the data frame and append
	newLine = data.frame(Species = species, Library = library, LINE = lines, LINE_perc = lines_perc, SINE = sines, SINE_perc = sines_perc, LTR = ltr, LTR_perc = ltr_perc, DNA = dna, DNA_perc = dna_perc, Unknown = unknown, Unknown_perc = unknown_perc, Total = total, Total_perc = total_perc)
	table = rbind(table, newLine)

}

boo = table$Species == "HD"
table$Species[boo] = "Hypsibius dujardini"
table$Species[!boo] = "Ramazzottius varieornatus"

write.table(x = table, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)