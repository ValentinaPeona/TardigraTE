# Practical 3: manual curation of consensus sequences

## Introduction
In the first practical you run RepeatModeler2 and obtained a library of raw consensus sequences for RV. We call them "raw" because most of the times, the precise sequence and termini of the repeats are not well definied especially when the repeat is very long or the genome assembly used is fragmented. To get more precise consensus sequences, the raw consensi need to be manually curated by following these steps:

1. create an alignment for each consensus sequence and its 20 best hits in the genome
2. generate a new consensus sequence (e.g., with Advanced Consensus Maker)
3. find the termini
4. curate the gaps, CpG sites and uncertain base calling
5. find the TSDs (if any)
6. document the curation and classification process (e.g., by using the Excel template given in this course)

Here you can find some code to get the alignments to curate for each consensus sequence. Briefly, we give you a Perl script `RMDL_curation_pipeline.pl` that uses Blast and MAFFT to create and filter the alignments.
The script needs only a fasta file with the consensi, the fasta file of the original genome assembly and a BLAST database for the genome assembly. At the end you'll find a folder called `final` with your alignments to curate.

---

# Part A

## Practical A.1: Run the alignment pipeline

Prepare your working directory

```bash
cd ~/TEcourse/Intermediate/
mkdir RV_rm1.0_curation && cd RV_rm1.0_curation
```

Collect the input files

```bash
# link the genome assembly
ln -s ~/TEcourse/Data/RV.fasta .

# link your file of consensi
GROUP=g01
ln -s ~/TEcourse/Intermediate/RMSK/RV_rm1.0_${GROUP}.lib .
```

**Important**: Before running the pipeline you need to make sure that no `/` characters are present in the headers of the repeat library.

```bash
awk '{print$1;}' RV_rm1.0_${GROUP}.lib | sed 's/\//_/g' > RV_rm1.0_${GROUP}.clean.lib
```

Make the BLAST database for the genome assembly

```bash
# make Blast database before running the Perl script
makeblastdb -in RV.fasta -out RV.fasta -dbtype nucl -parse_seqids
```

The code below may take too long but you can try to run it with a smaller fasta of consensus sequences. Then copy the all the alignments for your group from the `Share/Output/RV_rm1.0_curation/` folder

Subset your repeat library

```bash
perl -pe '/^>/ ? print "\n" : chomp' RV_rm1.0_${GROUP}.clean.lib | tail -n +2 | head -n 10 > RV_rm1.0_${GROUP}.sub.lib
```

```bash
REF=RV.fasta
CONS=RV_rm1.0_${GROUP}.sub.lib

# run the RModeler_pipeline.pl script
perl ~/TEcourse/Code/RMDL_curation_pipeline.pl $REF $REF $CONS

# remove temporary files
rm *emp.out

# create output folder
mkdir final

# reshape files
cd aligned; for i in $(ls *.fa); do name=`ls $i | cut -f1 -d "."`; cat $i | perl -ne 'chomp;s/>\s+/>/;if(/>(\S+)/){$id{$1}++;$id2=$1;}if($id{$id2}==1){print "$_\n"}' >../final/$name.fa; done; cd ../

# remove alignents that are too "gappy"
cd final; for i in $(ls *.fa); do name=`ls $i | cut -f1 -d "."`; t_coffee -other_pg seq_reformat -in $i -action +rm_gap 95 >$name.gaps95.fa; done; cd ../
```

The code generated 3 different folders

- `blast`: contains the BLAST alignments of your repeats in the genome assembly
- `aligned`: contains the alignments (done by MAFFT) between the occurrences of the repeats (found by BLAST) to the raw consensus sequences
- `final`: contains the filtered alignments from aligned.

In the `final` folder you will find two types of fasta alignments for each consensus sequence: one just called `.fa` and one `.gaps95.fa` in which sequences that regions with too many gaps have been removed. It is important to keep them both in case you want to retrace what happened to your alignments.

You can find all the alignments for your groups in `~/Share/Output/RV_rm1.0_curation/` and on Dropbox https://www.dropbox.com/sh/gqj8nrfvdav5rhh/AADYHJnpZa5GocWCKUFn7yG9a?dl=0

**Important**: Please try to have a Dropbox folder on your computer that synchronises with the shared Dropbox folder of the course. In this way it is easier to avoid file misplacements and redundancies. If you don't manage to have a synchronising Dropbox folder, please be careful in uploading your files into the right subfolders!!

On Dropbox in the folder physaliaTEcourse2021 > Collaborative_project you can find subfolders associated to your groups. Inside there are the alignments ready for you.

**Important**: Within your group subfolder, create the following directories: 

```
# Example for group 1

Group1_curation
	|___ 01_complete 	# put the curated alignments here
	|___ 02_incomplete 	# put the alignements that needs further curation here
	|___ 03_check 		# put any alignment you cannot curate completely
	|___ 04_what_is_this 	# put the alignments that you don't know how to curate
	|___ 05_not_curatable 	# put the alignments that are not curatable

```


**Important**: Unzip the compressed file you find in your respective group folder on Dropbox. You'll find two additional folders `Raw` and `Alignments`. Curate the alignments within `Alignments` with the extension **.gaps95.fa**

**Important**: DO NOT CURATE THE ALIGNMENTS IN `Raw`

Give a look at the `.gaps95.fa` alignments by opening them with BioEdit, AliView (or your favorite software). Can you find the TE within the chaos?

Good luck with the curation :P

Dropbox link: https://www.dropbox.com/sh/urt48s9n2166rwn/AADiUC3_YY-n2NOce3AO351Ya?dl=0