# -*- snakemake -*-

# Pipeline to run the analysis for the TardigraTE 

#############################################################################
# ===========================================================================
# Valentina Peona
# Alexander Suh Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 08 December 2022
# ===========================================================================
#############################################################################

# Usage: snakemake --snakefile TardigraTE_1.3.smk --profile smk_profile/ --cluster-status smk_profile/status-scontrol.sh

# --- Importing Configuration Files --- #

# configfile: "TardigraTE_config.yaml"

# -------------------------------------------------

# reference = config["Reference"]
# species = [i.split('/', 1)[1].split(".fasta")[0] for i in reference]
# TElibraries_path = config["TElibraries_path"]
# TElibraries = config["TElibraries"]
# lib = config["LibPrefix"]
# library_dict = dict((zip(lib, TElibraries)))
# ref_dict = dict((zip(species, reference)))

reference = ["Data/HD.fasta", "Data/RV.fasta"]
assemblies = ["HD.fasta", "RV.fasta"]
species = [i.split('/', 1)[1].split(".fasta")[0] for i in reference]
# All the TE libraries are assumed to be in this folder: Data/Libraries
TElibraries = ["Repbase.lib", "Uncurated.lib", "Curated.lib"]
lib = [i.split('.lib', 1)[0] for i in TElibraries]
library_dict = dict((zip(lib, TElibraries)))
ref_dict = dict((zip(species, reference)))
assembly_dict = dict((zip(species, assemblies)))
# These are the libraries specific for the single species
# in the same folder (Results/Libraries) also "HD_uncurated_single.lib" and "RV_uncurated_single.lib"
# MUST be present
lib_curated_single = ["HD_curated_single.lib", "RV_curated_single.lib"]
lib_uncurated_single = ["HD_uncurated_single.lib", "RV_uncurated_single.lib"]
mask_dict = dict(zip(['RV', 'HD'], lib_curated_single))

# -------------------------------------------------

localrules: all, TableOverviewConsensi, TablePercentage, TablePutativeActive, FigureLengthComparison

# -------------------------------------------------

rule all:
	input:
		OUT = expand("Intermediate/RMSK/{species}_{lib}.out", species = species, lib = lib),
		TBL = expand("Intermediate/RMSK/{species}_{lib}.tbl", species = species, lib = lib),
		ALIGN = expand("Intermediate/RMSK/{species}_{lib}.align", species = species, lib = lib),
		DIV = expand("Intermediate/RMSK/{species}_{lib}.divsum", species = species, lib = lib),
		TABLE1 = "Results/Tables/TableOverviewConsensi.tsv",
		PDF = expand("Results/Figures/{species}_{lib}.pdf", species = species, lib = lib),
		TABLE2 = "Results/Tables/TablePercentage.tsv",
		FIGURE2 = "Results/Figures/FigureLengthComparison.pdf",
		TABLE3 = "Results/Tables/TablePutativeActive.tsv"

# -------------------------------------------------

rule RepeatMasker:
	input:
		PLACEHOLDER = reference
	output:
		OUT = "Intermediate/RMSK/{species}_{lib}.out",
		TBL = "Intermediate/RMSK/{species}_{lib}.tbl",
		ALIGN = "Intermediate/RMSK/{species}_{lib}.align"
	params:
		LIBS = lambda wcs: library_dict[wcs.lib],
		SPECIES = species,
		REF = lambda wcs: ref_dict[wcs.species],
		ASSEMBLY = lambda wcs: assembly_dict[wcs.species]
	resources:
		threads = 16,
		mem_mb = 96000
	shell:
		"""
		module load bioinfo-tools
		module load RepeatMasker/4.1.0
		
		echo {input.PLACEHOLDER}
		
		DIR=`echo {output.OUT} | sed 's/.out//' | sed 's%Intermediate/RMSK/%%'`
		RepeatMasker -pa 16 -a -xsmall -gccalc -excln -s -dir Intermediate/RMSK/$DIR -lib Results/Libraries/{params.LIBS} {params.REF}

		mv Intermediate/RMSK/$DIR/{params.ASSEMBLY}.out {output.OUT}
		mv Intermediate/RMSK/$DIR/{params.ASSEMBLY}.tbl {output.TBL}
		mv Intermediate/RMSK/$DIR/{params.ASSEMBLY}.align {output.ALIGN}
		"""

# -------------------------------------------------

rule TableOverviewConsensi:
	input:
		LIBS = expand("Results/Libraries/{library}", library = TElibraries)
	output:
		TABLE1 = "Results/Tables/TableOverviewConsensi.tsv"
	shell:
		"""
		module load R_packages/4.1.1

		Rscript --vanilla Code/TableOverviewConsensi.R {output.TABLE1} {input.LIBS}
		"""

# -------------------------------------------------

rule TablePercentage:
	input:
		TBL = expand("Intermediate/RMSK/{species}_{lib}.tbl", species = species, lib = lib)
	output:
		TABLE2 = "Results/Tables/TablePercentage.tsv"
	shell:
		"""
		module load R_packages/4.1.1

		Rscript --vanilla Code/TablePercentage.R {output.TABLE2} {input.TBL}
		"""

# -------------------------------------------------

rule CalcDivergence:
	input:
		ALIGN = "Intermediate/RMSK/{species}_{lib}.align"
	output:
		DIV = "Intermediate/RMSK/{species}_{lib}.divsum"
	shell:
		"""
		module load bioinfo-tools
		module load RepeatMasker/4.1.0

		calcDivergenceFromAlign.pl -s {output.DIV} -a {input.ALIGN}
		"""

# -------------------------------------------------

rule CreateRepeatLandscape:
	input:
		OUT = "Intermediate/RMSK/{species}_{lib}.out"
	output:
		PDF = "Results/Figures/{species}_{lib}.pdf"
	resources:
		threads = 1,
		mem_mb = 6000,
		time = "00:15:00"
	shell:
		"""
		module load R_packages/4.1.1
		
		# determine assembly siye for this run
		TBL=`echo {input.OUT} | sed 's/.out/.tbl/'`
		SIZE=`grep 'total length' $TBL | tr -s ' ' | cut -f5 -d ' ' | sed 's/(//'`
		
		Rscript --vanilla Code/CreateRepeatLandscape.R {input.OUT} {output.PDF} $SIZE
		"""

# -------------------------------------------------

rule TablePutativeActive:
	input:
		OUT = expand("Intermediate/RMSK/{species}_{lib}.out", species = species, lib = lib)
	output:
		TABLE3 = "Results/Tables/TablePutativeActive.tsv"
	shell:
		"""
		module load R_packages/4.1.1

		Rscript --vanilla Code/TablePutativeActive.R {output.TABLE3} {input.OUT}
		"""

# -------------------------------------------------

rule FigureLengthComparison:
	input:
		CURATED = expand("Results/Libraries/{library}", library = lib_curated_single),
		UNCURATED = expand("Results/Libraries/{library}", library = lib_uncurated_single)
	output:
		FIGURE2 = "Results/Figures/FigureLengthComparison.pdf"
	shell:
		"""
		module load R_packages/4.1.1

		Rscript --vanilla Code/FigureLengthComparison.R {output.FIGURE2} {input.CURATED} {input.UNCURATED}
		"""

# -------------------------------------------------
