# functional and ecological determinants of evolutionary dynamics in crAss-like phages

contains bash scripts for 
* running annotation pipeline 
* main pipeline: determining selective pressure and micro diversity for genes
** metagenomic read mapping 
** data analysis

NB: when in doubt, consider using subversion commit history with commit messages part of the documentation

## examples ##
Accompanying Python scripts are here:
https://github.com/resharp/phages


## annotation pipeline ##

NB: scripts in annotation pipeline are run individually for maximum control and flexibility

-	Guerin_to_ref.sh
	-	split Guerin.fna in different fasta files in different subdirectories for all reference genomes
	-	predict genes and proteins on reference genomes with prodigal
-	Yutin_to_hmm.sh
	-	download alignments
	-	build HMM profiles for gene families
	-	create protein file for genus 1 from GenBank file
	-	search all reference genomes with HMM profiles, outputs:
		-	ref_genes.tsv
		-	hmm_hits.tsv
	-	run a blastp search of all proteins against all proteins
-	run_ref_against_pvog_hmm.sh
	-	run HMM search on all genes against HMM profiles in pVOG database (Grazziotin 2017, date accessed 14/1/2020)
-	folder hmm_hits
	-	hmm_hits.tsv
		-	results of hmm hits on correct gene annotation for genomes with different genetic code
	-	ref_genes.csv
		-	result of new gene predictions with correct translation table for two reference genomes
	-	pvogs_annotations.tsv
		-	pvogs annotations (by Nikos Pappas) from pVOG database (Grazziotin 2017, date accessed 14/1/2020)
	-	all_refs_pvog_table.txt	
		-	all pvog hits against ref genome genes
	- 	pairwise_blastout_filtered.txt
-	folder ref_seqs
	-	ids_ref_genomes.txt
		-	ref genomes that are used by run_pipeline_on_project.sh
		-	contains mapping from name reference genome to genus
	-	crassphage.gbk
		- original gbk file and results from create_proteins_from_gbk added.
	-	crassphage_refseq.fasta
		reference sequence
	-	crassphage_refseq_codingregions.txt
		-	coding regions for crassphage_refseq (used as input for DiversiTools)
	-	crassphage_refseq_gene_list.txt
		-	annotation outline output
	-	IDEM for all other reference sequences

## main pipeline: determining selective pressure and micro diversity for genes ##

## metagenomic read mapping ##
Main pipeline (ran for 3 weeks on 800 files) preparing for data analysis.

-	run_adapter_removal.sh
	-	Adapter removal and quality filtering of the paired-end fastq files was done with AdapterRemoval v2
-	run_pipeline_on_project.sh
	-	run_all_samples
		-	main pipeline, ran for three weeks for 800 samples!
		-	executes metagenomic read mapping, DiversiTools and CalcDiversiMeasures.py

temporary
-	run_pipeline_untrimmed
	-	used to determine if results differ for trimmed or untrimmed files

## data analysis pipeline ##

Wrappers/examples for MakeGenePlots.py, CalcDiversiMeasures.py and CalcCodonMeasures.py
Runs in the order of minutes!

-	data_analysis.sh
	-	run_all_codon_measures
	-	run_all_gene_plots
        -	analysis for one reference genome (individual gene analysis)
	-	run_family_gene_plots
		-	analysis for gene families based on multiple ref genomes
	-	sample_stats
		-	collect sample/genome mapping statistics produced by samtools idxstats
		-	output: sample_stats.txt (input file for MakeSamplePlots.py

## obsolete code used for pilot
Samples were selected from a local database MGXDB
-	run_crassphage_samples.sh
-	run_diversitools.sh


# old scripts of the phage clustering pipeline internship Bas Dutilh



## extracting phages ##
Phages were extracted using VirSorter from a selection of 34590 bacterial genomes from PATRIC.
https://www.patricbrc.org/

Use the script run_parallel_vs.sh to start 24 jobs of run_vs_1.sh in parallel distributed over two servers.
run_vs_1.sh uses the VirSorter Perl script and is set up in a Conda environment provided by VirSorter.

See
https://github.com/simroux/VirSorter

job_puller, ana_iterations and show_virsorter_summary are old script files, for running on concatenated PATRIC genomes.


## gene annotation ##

* We use prodigal for prediction of genes
* We do a pairwise blastp search of all genes against all genes with evalue cut-off 1e-3
* We filter out all pairwise matches with a > 75% query and target coverage
* we did _not_ filter on bitscore (compare vCONTact bitscore > 50)


## protein clustering ##
* We use mcl for clustering the genes based on the filtered pairwise blastp matches without weight on the edges
* mcl with inflation factor 2.0
* we evaluate the protein clusters by building a hmm profile for every gene cluster and see how well it separates the proteins (todo: specify)
* mafft for multiple alignment of proteins
* hmmbuild for building hmm profile
* hmmsearch

## clustering phages ##
* we cluster on shared gene content

