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


## gene annotation ##

* We will use prodigal for prediction of genes


## clustering phages ##







