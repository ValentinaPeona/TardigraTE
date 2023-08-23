# Practical 1: De-novo repeat prediction and annotation with RepeatModeler2 and RepeatMasker

*Valentina Peona*

*valentina.peona@ebc.uu.se*

## Introduction
Masking a genome of a new species using repeat libraries from other species can be insufficient. It is particularly important when the genome of interest contains species-specific repetitive elements never described before and these elements would remain largely unmasked or partially and/or incorrectly annotated. Furthermore, curated repeat libraries from sister species can actually reciprocally help the annotation of the two genomes (Boman et al. 2018).

To get a more complete annotation we need to run a de-novo repeat annotation on our genome of interest to get a new repeat library made of raw consensus sequences. These raw consensus sequences are going to need a good round of manual curation and after that the library will be ready to annotate our genome of interest. At the end of our practicals we'll compare the genome annotation made with raw consensus sequences and with curated consensus sequences; the difference is going to be quite remarkable!

There are a few tools available for the de-novo characterisation of repeats (e.g., CARP, REPET, RepeatExplorer2 and other all listed in https://tehub.org) but for time reasons we will use only RepeatModeler2. As target genome, we chose to use the recently assembled genome of the tardigrade *Ramazzottius varieornatus* (abbreviated throughtout the practicals as "RV"; Yoshida et al. 2017). Once the RV repeat library is complete we will run RepeatMasker with the new library on the RV genome.

In the previous edition of this course, your colleagues curated the repeat library of another tardigrade *Hypsibius dujardini* (abbreviated as "HD"; Yoshida et al. 2017) and this allows us to also compare the repeat libraries of two different species and gives us the opportunity to learn how to merge separate repeat libraries in an efficient and non-redundant way.

### Note
Given the time we have for the course, all the commands have been already run so you can find all the outputs in the folder `~/Share/Output/`. We suggest you to follow the tutorial and try the commands out but full genome runs of RepeatModeler and RepeatMasker would take way too much time therefore kill the jobs once they started correctly.

In this Practical we focus on:
1. running RepeatModeler2 on RV genome
2. inspecting the output of RepeatModeler2
3. running RepeatMasker with the new library
4. inspecting the output of RepeatMasker

### Funny note
The genus Ramazzottius was not named after the Italian singer Eros Ramazzotti but after the Italian naturalist Giuseppe Ramazzotti. He was part of the family that produces the Amaro Ramazzotti liquors but decided to dedicate his life to science instead of to drinks and money ò.ò if and when you get a chance to drink this liquor, remember to cheer for our beloved tardigrades \^^/

---
# Part A

## Practical A.1: Install RepeatModeler2 (optional)

On the Amazon server RepeatModeler (Flynn et al. 2020) and RepeatMasker (Smit, Hubley, Green 2010) have been installed for you but here some instructions for installing it on your own server/computer.

**Follow these instructions ONLY if you need to install RepeatModeler2 on your computer!**

The new version of RepeatModeler has been released in 2020 and published in PNAS https://www.pnas.org/content/117/17/9451.short (Flynn et al. 2020) and can be downloaded from the GitHub page https://github.com/Dfam-consortium/RepeatModeler

There are several ways to install RepeatModeler2, for example through Conda or by directly downloading a Singularity image of the software which contains all the dependencies needed.

**Conda installation**
```bash
conda create -n RepeatModeler2
conda activate RepeatModeler2
conda install -c bioconda repeatmodeler
```

Once the environment has been created correctly, you can run RepeatModeler2 just by calling it on the command line.

There are several other tools you can install alongside RepeatModeler to run the LTRStruct option (to improve LTR prediction), please follow the instructions at https://github.com/Dfam-consortium/RepeatModeler

**Singularity installation**

A Singularity and Docker image can be downloaded from https://github.com/Dfam-consortium/TETools

Example of how to install it and run

```bash
cd /path/to/your/folder/
curl -sSLO https://github.com/Dfam-consortium/TETools/raw/master/dfam-tetools.sh
chmod +x dfam-tetools.sh
./dfam-tetools.sh
```

Inside the container, the included tools are now available:

```bash
BuildDatabase -name genome1 genome1.fa
RepeatModeler -database genome1 [-LTRStruct] [-pa 8]
```

If neither of these installation options work for you, please visit the Github pages: https://github.com/Dfam-consortium/RepeatModeler and https://github.com/Dfam-consortium/TETools.

---

## Practical A.2: Getting ready

**Important bash commands**

```bash
# change directory
cd

# get the name of the folder in which you are located
pwd

# list files in directory
ls
ll

# visualise a file
less <filename>
# exit less by typing q

# print the content of a file
cat <filename>

# print first lines of a file
head <filename>

# print last lines of a file
tail <filename>

# edit a file
nano <filename>

# create directory
mkdir <dir_name>

# copy file
cp /path/to/file/file.txt ./

# create symbolic link
ln -s /path/to/file/file.txt ./

# get help for how to run commands
man <command>

# if man does not work try with
<command> -h
```

**Create your working directory**

We need to collect all the input files and scripts necessary for this course and copy them to your own folders:

1. the *Ramazzottius varieornatus* genome `RV.fasta`
2. the curated repeat library for *Hypsibius dujardini* `HD.lib`
3. Perl script to rename the consensus sequences created by RepeatModeler2 `renameRMDLconsensi.pl`
4. Perl script to recalculate the divergence from consensus `calcDivergenceFromAlign.pl`
5. Perl script to create landscape plots `createRepeatLandscape.pl`
6. Perl script to generate alignments to manually inspect `RMDL_curation_pipeline.pl`
7. Some Perl and Python scripts to manipulate RepeatMasker output files `RM2Bed.py` `rmOutToGFF3.pl` `rmOut2Fasta.pl`

**Note** all the scripts (but `RMDL_curation_pipeline.pl` and `renameRMDLconsensi.pl`) are part of the RepeatMasker installation BUT here on this server it seems complicated to reach them, therefore we provided a local copy of them. They can be usually found in the folder `RepeatMasker/util/` on your own server/computer.

All the files and scripts are given in the folder `Share` in your own home directory. Please create your working directory and copy the files from `Share`.

Throughout the tutorial a folder system that looks like this will be used but you're free to use the system you prefer:

```
TEcourse (main working directory)
        |_________Data # contains the original input files
        |_________Intermediate # contains all the output files from the analysis
        |_________Results # contains all the final outputs
        |_________Code # contains code files and scripts
```

Create your own folders

```bash
# create folders
mkdir -p TEcourse/Data TEcourse/Code TEcourse/Intermediate TEcourse/Results

# copy reference genome and repeat library in the Data folder
cp ~/Share/Data/RV.fasta ./TEcourse/Data/
cp ~/Share/Data/HD.lib ./TEcourse/Data/

# copy scripts in the Code folder
cp ~/Share/Code/* ~/TEcourse/Code/
```

A copy of all the files can also be found on Dropbox in the folder AWS_server_files: 
https://www.dropbox.com/sh/xn6snum1853jf32/AAAdCyMctjz-HplGeb9I9tC1a?dl=0

---

## Practical A.3: Run RepeatModeler2

On the Amazon server we installed RepeatModeler2 through conda and you should be able to run it by simply calling the command RepeatModeler from the command line.

Once you have your files copied we can start giving a look at RepeatModeler2 and its help page:

```bash
# create output folder for RepeatModeler in Intermediate
conda activate repeatmodeler

cd ~/TEcourse/Intermediate
mkdir RMDL
cd RMDL

RepeatModeler -h
```

Inspect the options, do you see anything useful/interesting? Which options would you use? What does `LTRStruct` mean/do? (The installation of the dependencies for LTRStruct resulted to be quite hard on the Amazon server so we will not use it in the command below BUT it is much easier to install on your own computers!)


### Note
**Please run the first command and start the second but then kill it and copy the output we produced for you in your folder.**

Let's give it a try

```bash
# symlink the genome file in your RMDL folder
ln -s ../../Data/RV.fasta .

# STEP1: create a database for the genome
BuildDatabase -name RV RV.fasta

# STEP2: run RepeatModeler2 on the database you just created
RepeatModeler -database RV
```

**Kill the RepeatModeler command, we created the output files for you already** kill command: `CTRL + C`

As you can see there are two commands to run to get the final output from the tool (the library of raw consensus sequences). The first command creates a database of your genome assembly fasta file, basically it indexes the fasta to better access it during the de-novo characterisation similarly to what BLAST. Indeed, most of the RepeatModeler analysis consists of a large number of alignments.

After that the database has been correctly generated, the second command will start the true analysis.


### Note
When you're going to run RepeatModeler and RepeatMasker on your fasta files, pay attention to the presence of strange/special characters in your fasta headers, they may cause problems during the analysis. Also it is suggested to have rather short fasta headers or they will be cut and you may lose some important pieces of information: this is of particular importance when running RepeatMasker.

**Copy the output files in your own RMDL folder**
```bash
cp -r ../../../Share/Output/RMDL/* .
```

---

## Practical A.4: Inspect the RepeatModeler2 output

The output is rather straightforward to understand. In the output folder we can find four key files: `consensi.fasta`, `consensi.fa.classified` and two Stockholm files for the two previous files.

The `consensi.fasta` contains the raw consensus sequences named by RepeatModeler automatically.

The second and more important file is the `classified` version of the consensus sequences. One of the processes of the RepeatModeler analysis tries to preliminarly classify the sequences by aligning the consensus sequences to a database of known transposable element related proteins. This step is important to get a first overview of the type of transposable elements found in the genome but it can produce artifacts especially when species never analysed before are concerned. Misclassification is unfortunately typical and that's one of the reasons why we need to curate these consensus sequences in the next days!!

The last two files are Stockholm alignment files in which for each consensus sequence, seed alignments througout the genome are collected.

Give yourself a look at the files now :)
How many consensus sequences have been produced?
Which kind of repetitive elements do you find in the library?
Are the consensus sequences long, short?

---

# Part B

## Practical B.1: Install RepeatMasker (optional)


**Install RepeatMasker ONLY if you're working on your own server/computer**

For RepeatMasker there are no Singularity images available but we can always rely on Conda

```bash
# if you already have the RepeatModeler environment activated please deactivate it:
# conda deactive RepeatModeler
conda create RepeatMasker
conda activate RepeatMasker
conda install -c bioconda repeatmasker
```

and we're ready to go!

---

## Practical B.2: Compare RV and HD libraries with RepeatMasker

In this part of the tutorial we learn how to run RepeatMasker in groups of two or three people. **Be sure you are in the right breakout room and know your group name.**

Since running the entire new library on the genome assembly will take likely a couple of hours, we're going to run it using smaller files.
We're going to mask the entire HD curated repeat library (reference) with a subset of the RV raw repeat library (query). In this way, in a little time, we discover if there is homology between the RV and HD repeats. *The results from Practical B.2 will be the basis for practicing library merging in Practical 1 part D.*

In the folder `Share/Output/RMDL/Groups`, you can find your subset of RV consensus sequences named after the number of your groups. Copy your file into your working directory.

**Note: please use your own group number below!!**

```bash
# create RepeatMasker output folder in Intermediate
# go to Intermediate
cd ~/TEcourse/Intermediate
mkdir RMSK
cd RMSK

# example using group 1; change the GROUP variable according to your group name
# group 2: g02; group 3: g03; group 10: g10 etc etc
GROUP=g01
cp ~/Share/Output/RMDL/Groups/RV_raw_consensi_${GROUP}.fasta .

CODE=~/TEcourse/Code

# rename consensus sequences: replace the standard "rnd" suffix with a custom one
# ramVar = Ramazzottius varieornatus
perl $CODE/renameRMDLconsensi.pl RV_raw_consensi_${GROUP}.fasta ramVar RV_rm1.0_${GROUP}.lib

# link the curated repeat library for HD tardigrade HD.lib from Data
ln -s ../../Data/HD.lib .
```

Before running RepeatMasker, let's give a look at the different options.
List the options by typing:

```bash
# go back to the base conda env
conda deactivate
conda activate RepeatMasker

RepeatMasker
#or for having a more detailed help
RepeatMasker -help
```

Since the manual is quite long, you can temporarly save it in a file

```bash
RepeatMasker -help > RMSK_man
less RMSK_man
```

Some important options you should look at are

```
 # number of cores to use
 -pa(rallel) [number]

 # options for the speed of the repeat annotation
 -s 	slow *RECOMMENDED*
 -q 	quick
 -qq 	very quick

 # options for including or excluding some repeat categories
 -nolow 	do not mask low complexity
 -noint 	do not mask interspersed repeats - mask only simple/low complexity
 -norna 	do not mask RNA genes
 -alu 		mask only Alu (relevant for primates)
 -div 		mask repeats below a certain divergence threshold (e.g., useful when you want to know only about young repeats)

# option for using a specific repeat library!!
-lib 	specify the path to your custom library *RECOMMENDED*

# option to specify a species name
-species 	RepeatMasker will use consensus sequences specific to that species (e.g., if specify mouse, it will use sequences specific to mouse, sequences present in all Rodents and in all Mammals). -lib and -species are mutually exclusive *RECOMMENDED*

# option for dealing with contamination (bacterial insertion sequences)
-is_only
-is_clip
-no_is

# options for masking
-gccalc 	calculates the GC content of the sequence to better distinguish repeats from low complexity sequences *RECOMMENDED*

# options for the output
-dir 		you can specify an output folder *RECOMMENDED*
-a 		write alignment file (.align) *RECOMMENDED*
-xsmall		returns a fasta file with soft-masked repeats *RECOMMENDED*
-x 		returns a fasta file with hard-masked repeats
-html 		returns the RMSK output as an additional html file
-gff 		returns the RMSK output as an additional GFF file *IMPORTANT*
-e(xcln) 	outputs a table file with a summary of the repeat proportions *RECOMMENDED*
```

Which options would you use?

This is the command we're going to use.

```bash
# be sure you're using the right GROUP value
GROUP=g01
RepeatMasker -pa 1 -a -xsmall -gccalc -excln -s -dir ./ -lib HD.lib RV_rm1.0_${GROUP}.lib
```

---
# Part C

## Practical C.1: Inspect and use RepeatMasker output

In the output directory indicated above (`-dir`) you will find several output files:

#### .out
the main and most important output file of RepeatMasker. This files contains 15 columns describing in detail each repeat found in the genome.

Which is the column containing the divergence from consensus?
What do the columns 8 and 14 mean?

#### .align
the `.align` file contains all the alignments found by RepeatMasker and give some supplementary information about the alignment statistics. This file is important to keep for some downstream analysis below.

#### .tbl
this is a summary table (the one produced thanks to the option `-excln`) containg useful information about the portion of sequence masked and proportion of repeats.

#### .masked
this is the soft-masked version of the input fasta file

#### .cat.gz
this is the version of the `.align` file before the `ProcessRepeats` step of RepeatMasker

Now we explore some useful tools and commands to investigate and manipulate the RepeatMasker output files. Many of the scripts listed below are found in the `util` folder within the installation folder of RepeatMasker itself.

### Re-calculate the divergence adjusted for GC content

The presence of methylated cytosines at CpG sites (in animals) often causes the hypermutability of Cs into Ts (leading to TpG and CpA sites) and this mechanism can cause an overestimation of the percentage of divergence from the consensus sequence. Therefore, there is a Perl scripts in the utils of RepeatMasker that can re-calculate the divergence percentages taking into account the presence of CpG sites starting from the `.align` file.

**Note** the CpG re-calculation given here is heavily biased for animals and other organisms in which the context CpG is methylated so if you are using your own data and you know CpGs are not methylated you can skip this step. Also it should be possible to modify the script to consider other methylation contexts.

Here 

```bash
CODE=/home/ubuntu/miniconda3/envs/RepeatMasker/share/RepeatMasker/util/

cd ~/TEcourse/Intermediate/RMSK/
mkdir RV
cd RV

cp ~/Share/Output/RMSK/RV_rm1.0/RV.fasta.out .
cp ~/Share/Output/RMSK/RV_rm1.0/RV.fasta.align .

# Calculate the divergence with Kimura 2-parameter distances (excluding CpG sites)
perl $CODE/calcDivergenceFromAlign.pl -s RV.fasta.align.divsum -a RV.fasta.align_with_div RV.fasta.align
```

What is `-s`? What is `-a`?

Some more info about the CpG sites directly from the help of the `calcDivergenceFromAlign.pl` script:

Treat "CG" dinucleotide sites in the consensus sequence as follows:
Two transition mutations are counted as a single transition, one transition is counted as 1/10 of a standard transition, and transversions are counted normally (as the would outside of a CpG site).  This modification to the Kimura 2 parameter model accounts for the extremely high rate of mutations in at a CpG locus.

Look at either of the new files: Which is the column containing the new divergence estimate? Is it in the same format as the previous divergence? Did the divergence change? If it changed, did it increase or decrease?

### Make a repeat landscape plot

To make a nice and fast landscape plot we can easily use another RepeatMasker script that will take the output from the previous step as input.

```bash
CODE=/home/ubuntu/miniconda3/envs/RepeatMasker/share/RepeatMasker/util/

# Make repeat landscape from CpG-corrected RepeatMasker .align file and average divergence file:
perl $CODE/createRepeatLandscape.pl -div RV.fasta.align.divsum -g 55423593 > RV.fasta.align.divsum.html
```

What is `-g`?
The number given in the command line is specific for the RV genome, when running the command on another genome you must change that parameter according to the genome size! You can find the size of the fasta analysed in the `.tbl` file


Download locally the `.html` file and give it a look! Isn't it pretty? :D

### Convert an align file into .out (optional)

From the `calcDivergenceFromAlign.pl`, we obtained a new `.align` file with the CpG-corrected divergence but it would be nice to have it in the form of the `.out` file again so here's a one-liner awk command to do so:

```bash
awk '{if(index($12, "(") != 0){print $1,$2,$3,$4,$5,$6,$7,$7-($6-1),$8,"+",$9,$10,$11,$12,$13,$14}else if($9 ~ /C/){print $1,$2,$3,$4,$5,$6,$7,$7-($6-1),$8,$9,$10,$11,$12,$13,$14,$15}}' OFS="\t" RV.fasta.align_with_div | sed 's/#/\t/' > RV.fasta.align_with_div.noCpG.size.out
```

### Convert the .out file into a BED file (optional)

Sometimes it's very handy to have a BED file to use with BEDTools (e.g., intersect repeat annotation with other types of annotations; Quilan and Hall 2010).
We show you three ways to convert the `.out` into BED

1. Make 0-formatted bed file from 1-formatted RepeatMasker .out file (sequence, masking start, masking end) with awk:

```bash
awk '(NR>3){$6=$6-1; print $5,$6,$7}' OFS="\t" RV.fasta.out > RV.fasta.out.bed
```

2. Make bed6 file from RepeatMasker .out file (bed file with TE information):

```bash
awk '(NR>3){$6=$6-1; print $5,$6,$7,$10,$11,$9}' OFS="\t" RV.fasta.out > RV.fasta.out.bed6
```

3. Make a BED file using a RepeatMasker util script

```bash
$CODE/RM2Bed.py RV.fasta.out
```

This script has several options, give them a look!
It also works starting from a `.align` file as input!

### Convert the .out file into a GFF file (optional)

GFF files are as useful as BED files therefore let's see how to convert the .out file into a GFF using a RepeatMasker util script

```bash
$CODE/rmOutToGFF3.pl  RV.fasta.out > RV.fasta.out.gff
```

### Get sequences from the .out file (optional)

The script `rmOut2Fasta.pl` is an easy and fast way to get a fasta file from the .out file

```bash
$CODE/rmOut2Fasta.pl -fasta RV.fasta -out RV.fasta.out
```
It can also be done using a BED file (+ the original fasta file) and the function bedtools getfasta from BEDTools.

### Make a hard-masked (NNN) version of the soft-masked (lowercase) .masked file (optional)

Sometimes it is useful to have a hard masked version of the genome and we can do this with a simple perl command:

```bash
perl -e 'while(<>) { if ($_ =~ /^>.*/) { print $_; } else { $_ =~ tr/acgt/N/; print $_;}}' < RV.fasta.masked >RV.fasta.masked.hard
```

You can also tune the soft and hard masking of a genome by using a custom BED track of repeats and `bedtools maskfasta`

### Programs in RepeatMasker for downstream analyses:
There are many other scripts included in the RepeatMasker installation that are useful for downstream analysis. Feel free to explore them!

**Scripts in the main RepeatMasker installation folder:**
- DateRepeats
- DupMasker
- ProcessRepeats
- RepeatMasker
- RepeatProteinMask

**Scripts in RepeatMasker/util folder:**
- buildRMLibFromEMBL.pl
- buildSummary.pl
- calcDivergenceFromAlign.pl
- createRepeatLandscape.pl
- dupliconToSVG.pl
- getRepeatMaskerBatch.pl
- queryRepeatDatabase.pl
- queryTaxonomyDatabase.pl
- RM2Bed.py
- rmOut2Fasta.pl
- rmOutToGFF3.pl
- rmToUCSCTables.pl
- trfMask
- wublastToCrossmatch.pl

---
## References

Boman, J. et al. The genome of blue-capped cordon-bleu uncovers hidden diversity of LTR retrotransposons in zebra finch. Genes 10, 301 (2019).

Flynn, J. M. et al. RepeatModeler2 for automated genomic discovery of transposable element families. Proceedings of the National Academy of Sciences 117, 9451-9457 (2020).

Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841-842 (2010).

Smit, A., Hubley, R. & Green, P. RepeatMasker Open-3.3.0. http://www.repeatmasker.org.  (1996-2010).

Yoshida, Y. et al. Comparative genomics of the tardigrades Hypsibius dujardini and Ramazzottius varieornatus. PLoS Biol. 15, e2002266 (2017).
