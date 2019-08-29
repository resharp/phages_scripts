# phages_scripts
bash scripts for the phage clustering pipeline internship Bas Dutilh

Richard Gremmen
gremmen@resharp.nl


## extracting phages ##
Phages were extracted using VirSorter from a selection of 34590 bacterial genomes from PATRIC.
https://www.patricbrc.org/

Use the script run_parallel_vs.sh to start 24 jobs of run_vs_1.sh in parallel distributed over two servers.
run_vs_1.sh uses the VirSorter Perl script and is set up in a Conda environment provided by VirSorter.

See
https://github.com/simroux/VirSorter

job_puller, ana_iterations and show_virsorter_summary are old script files, for running on concatenated PATRIC genomes.


## gene annotation ##

* We will use prodigal for prediction of genes
* We do a pairwise blastp search of all genes against all genes with evalue cut-off 1e-3
* We filter out all pairwise matches with a > 75% query and target coverage



## protein clustering ##
* We use mcl for clustering the genes based on the filtered pairwise blastp matches without weight on the edges
* mcl with inflation factor 2.0
* we evaluate the protein clusters by building a hmm profile for every gene cluster and see how well it separates the proteins (todo: specify)
** mafft for multiple alignment of proteins
** hmmbuild for building hmm profile
** hmmsearch

## clustering phages ##
* we will cluster on shared gene content
