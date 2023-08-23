# Practical 1: De-novo repeat prediction and annotation with RepeatModeler2 and RepeatMasker

*Valentina Peona*

*valentina.peona@ebc.uu.se*


## Introduction

In this practical you're going to reclassifly (if possible) the RV consensus sequences on the basis of the homology with HD repeats using the RepeatMasker output you produced yesterday at the end of Part B. This reclassification is part of the "TE field work" in which we try to assign a sample to its right taxonomical family (similar to genus) and subfamily (similar to species). The search for the right and complete taxonomy will actually end when all the sequences are curated, this is just a first step.

Usually, to define that two consensus sequences belong to the same family, the 80-80-80 rule is followed: "the “80-80-80” rule (Wicker et al. 2007) has been proposed as a means of identifying copies from the same TE family: two TE copies may be considered to belong to the same family if they are aligned, with 80% identity, over at least 80 bp and 80% of their respective lengths" (Flutre et al. 2011).

In addition to that, we can try to identify repeats that belong to the same subfamily with more stringent similarity parameters following the “95-80-98” rule. "This implies that a consensus sequence is *removed* if it is included within another consensus sequence over 98% of its length, with an identity level exceeding 95%" (Flutre et al. 2011).

---

# Part D


In this tutorial you're going to inspect the RepeatMasker output files that you produced yesterday at the end of Practical 1 Part B. You're going to check if any RV repeats is similar to HD repeats following the “95-80-98” (same subfamily; Flutre et al. 2011) or the “80-80-80” (same family; Wicker et al. 2007) rule.

Go to the folder in which you have the RepeatMasker output file for RV repeats vs HD repeats and inspect your file.

Initially inspect the file manually just to get a feeling of how similar RV repeats are to HD repeats then you can test the suggested code below.

```bash
cd ~/TEcourse/Intermediate/RMSK/

# change the group name!!
GROUP=g01
less RV_rm1.0_${GROUP}.lib.out
```

- How long are the homologous regions? Are they longer than 80 bps?
- What's the sequence similarity? (column 2)
- How many are similar to HD repeats?

It is essential for this part of the course and analysis that you record the results of Part D in the Excel template we made for you. Download Excel template from Dropbox: https://www.dropbox.com/s/8blfl16ztrylhc1/Classification_and_curation_factsheet_template.xlsx?dl=0

**TO DO!!!**: Paste all the names of your RV sequences into this template under the column O (rm1.0_name) and record similarity to the HD repeats in column Q. If there is no similarity to the HD repeats, also record this in column Q.

Example of how to print sequences names ready to paste into the Excel file:

```bash
grep '>' RV_rm1.0_${GROUP}.lib | cut -c2- > list_sequences
cat list_sequences
```

For the cases in which the homology requirements for family and/or subfamily are fullfilled, record the suggested reclassification in column R.

Put the resulting Excel file into the Dropbox folder: https://www.dropbox.com/sh/8d8008y0tudmtd9/AACzJd2kgJO335DYFDsdMR-ja?dl=0 **under your group subfolder**

**We're going to keep working on these files in the next days during the manual curation process**

**Optional**: Paste all sequences in CENSOR (https://www.girinst.org/censor/index.php) and see if there is similarity between RV consensus sequences and known repeats deposited in Repbase?

Before pasting the sequences please replace the "#" in the headers with another symbol because CENSOR is not happy about the "#" character. You can do it like this:

```bash
# replace # with double underscores
sed 's/#/__/' RV_rm1.0_${GROUP}.lib > RV_rm1.0_${GROUP}.lib.nohash
```

It is always better to keep all the file versions so do not overwrite your files.

**Optional code**

Here we propose a couple of awk commands to test the 80-80-80 rule.

The first awk command converts the .out file into a tab-separated file so it's easier to read and parse. The second command filters for the percentage of query covered by the alignment with HD repeats (higher than 80%). The last command filters for identity (higher than 80%)

```bash
cd ~/TEcourse/Intermediate/RMSK/

# change the group name!!
GROUP=g01

# convert .out file into a tab-separated file maintaining all the columns
awk 'NR >3 {print $0}' RV_rm1.0_${GROUP}.lib.out | tr -s ' ' | sed 's/^ //' | sed 's/[*]//' | tr -s ' ' '\t' > RV_rm1.0_${GROUP}.lib.tsv

# Filter for percentage of the query covered and with a cutoff of 80 bps
cat RV_rm1.0_${GROUP}.lib.tsv | sed 's/[()]//g' | awk '{print $0, $7-$6+1, $7+$8, ($7-$6+1)/($7+$8)}' | awk '$18 >= 0.8 {print $0}' > RV_rm1.0_${GROUP}.lib.tsv.len80

#Filter for identity above 80% (i.e., Kimura 2-parameter genetic distance of 0.2 and lower)
# RepeatMasker default output file but consider to use the CpG corrected file and change the column (after the $ symbol) for the divergence accordingly
awk '$2 <= 20' RV_rm1.0_${GROUP}.lib.tsv.len80 > RV_rm1.0_${GROUP}.lib.tsv.len80.k2p0.2
```

Don't worry if you don't find any homology at all, it's biology! :P Keep reading and you will find a way to train your re-classification skills on other files. 

**Note**: The commands above are for didactic purpose so all the steps are clearly stated. If you're fluent in awk/bash etc you're free to compact them.

Inspect the output files, is there any repeat that fulfils the "80-80-80" rule?

**Optional**: tweak the commands above to look for repeats fulfilling the “95-80-98” rule

As you may have noticed, the filtering commands depend on the masking direction (which library is the query and which one is the reference) and the format of the file (.out file vs BED, GFF etc) so please pay attention to the number of the columns when analysing other files. When running the filters also pay attention to which library is the curated one or to which consensus sequences are the longest etc (see slides of today's lecture).

This approach works quite well when we have a pair of libraries to compare and merge but how would you prioritise the library merging when dealing with multiple libraries? Which consensus sequences to keep/remove? How to reclassify them?

If you want to test this filtering/reclassification process on a raw repeat library from a species that has a closely related species with a curated repeat library (e.g., bird) you can give a look at the RepeatMasker output from 
raw Anna's hummingbird (*Calypte anna*) RMDL library vs the curated repeat library of several other birds.

Find the .out file for hummingbird vs paradise crow in: `Share/Output/RMSK/hummingbird/`
No need to record anything! It's just for training :)

```bash
cd ~/TEcourse/Intermediate/RMSK/

cp -r ~/Share/Output/RMSK/hummingbird/ .

cd hummingbird

# run the awk commands and look at the output
```

---

## Optional task 1: create and explore landscape plots from various libraries

Collect the RepeatMasker output files produced for you from the `Share/Output/RMSK` folder.

In `Share/Output/RMSK` folder you can find:
1. `RV_rm1.0`: RV genome masked with the RV raw consensus sequences (already analysed yesterda)
2. `Arthropoda`: RV genome masked with Repbase Arthropoda repeat library (collection of repeats known to be found in Arthropoda species)
3. `HDlib`: RV genome masked with the HD curated repeat library

Copy the files into you working directory and run `calcDivergenceFromAlign.pl` and `createRepeatLandscape.pl` as we did yesterday then compare the landscapes in the different `.html` files.

Example for `HDlib`:

```bash
cd ~/TEcourse/Intermediate/RMSK/

cp -r ~/Share/Output/RMSK/HDlib/ .

cd HDlib

CODE=/home/ubuntu/miniconda3/envs/RepeatMasker/share/RepeatMasker/util/
# Calculate the divergence with Kimura 2-parameter distances (excluding CpG sites)
perl $CODE/calcDivergenceFromAlign.pl -s RV.fasta.align.divsum -a RV.fasta.align_with_div RV.fasta.align

# Make repeat landscape from CpG-corrected RepeatMasker .align file and average divergence file:
perl $CODE/createRepeatLandscape.pl -div RV.fasta.align.divsum -g 55423593 > RV.fasta.align.divsum.html
```

**Note**: When running the commands for the different files, pay attention to the file names!

What differences can you find? Are there some categories that disappear/appear only in one landscape? What about the proportion of repeats at low divergence? What can you conclude?

---

## Optional task 2: run RepeatMasker on raw reads

You can run RepeatMasker directly on raw reads. It is interesting to compare the results on raw reads with the results obtained with the assembly.

In `Share/Output/Data` you can find two files of sampled raw reads of RV, one with a million reads and the other with half a million. Run RepeatMasker on them and create the landscape plots.

If for time reasons you don't manage to fully run RepeatMasker (it may take too long), you can find the .out files in 
`Share/Output/RMSK` subfolders `reads1M` and `reads500K`.

Example for reads1M:

```bash
cd ~/TEcourse/Intermediate/RMSK/

mkdir -p reads1M && cd reads1M

cp ~/Share/Data/reads1M/RV_1M_sampled.fasta .

conda activate RepeatMasker

ln -s ~/TEcourse/Intermediate/RMDL/RV_rm1.0.lib .

RepeatMasker -pa 1 -a -xsmall -gccalc -excln -s -dir ./ -lib RV_rm1.0.lib RV_1M_sampled.fasta
```

Create the landscape plot(s) and compare with the plot created using the assembly (follow the code of yesterday or the code above in Optional task 1). Which value should `-g` have for the calcDivergenceFromAlign.pl script? How do the landscapes differ?

Note that the reads have been directly sampled from the library without filtering steps. Which filtering steps would have you done? How would filtering affect the landscape plot? Would filtering steps affect the landscape at all?

---
# References

Flutre T, Duprat E, Feuillet C, Quesneville H (2011) Considering Transposable Element Diversification in De Novo Annotation Approaches. PLOS ONE 6(1): e16526. https://doi.org/10.1371/journal.pone.0016526

Wicker, T., Sabot, F., Hua-Van, A., Bennetzen, J. L., Capy, P., Chalhoub, B., ... & Schulman, A. H. (2007). A unified classification system for eukaryotic transposable elements. Nature Reviews Genetics, 8(12), 973-982. 