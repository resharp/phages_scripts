##
##
##
## example:
##
## out_dir=/hosts/linuxhome/mutant5/tmp/richard/virsorter_output/vs_all_1000
## show_virsorter_summary $out_dir >> your_output_file.txt
##
## ##OR do:
## show_vs $out_dir >> your_output_file.txt

function show_iteration_summaries {

	out_dir=$1

	echo "-----------------------------------"
	echo "Viral sequences found in iterations"
	echo "-----------------------------------"
	grep VIRSo -c $out_dir/r_*/global_signal_backup.csv

	echo "-------------------------------------------"
	echo "Growth of phage predictions with iterations"
	echo "-------------------------------------------"
	nr_iterations=$(grep -c VIRSo -c $out_dir/r_*/global_signal_backup.csv | wc -l)
	i=0
	while [ $i -lt $((nr_iterations-1)) ]
		do i=$((i+1))
		echo "$out_dir/r_$i/global_signal_backup.csv"
		grep -v "##" $out_dir/r_$i/global_signal_backup.csv | cut -d "," -f5 | sort | uniq -c
	done

	echo "$out_dir/VIRSorter_global-phage-signal.csv"

	grep -v "##" $out_dir/VIRSorter_global-phage-signal.csv | cut -d "," -f5 | sort | uniq -c
}

function determine_shared_genes {

	out_dir=$1

	echo "---------------------------------------------------------------"
	echo "Preparing the gene presence/absence matrix: determine pangenome"
	rm -f /dev/shm/phages_that_contain_phage_genes.txt
	rm -f /dev/shm/shared_phage_genes.txt
	rm -f /dev/shm/phages_plus_phage_genes.txt
	rm -f /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt
	rm -f /dev/shm/phages_with_shared_pangenome.txt

	##phages that contain phage genes
	join -1 1 -2 1 <(cat $out_dir/VIRSorter_global-phage-signal.csv | grep -v "^#" | cut -d ',' -f1| sort )  <(cat $out_dir/Metric_files/VIRSorter_affi-contigs.tab | grep -v '>' | tr '|' '\t' | awk '{print $1,"\t",$6}' | awk '{if($2!="-" && $2!="")print $0}' | sed -e 's/\-gene_[0-9]*//g' | cut -f1 | sort | uniq) \
		>> /dev/shm/phages_that_contain_phage_genes.txt

	##make list of shared phage genes
	join -1 1 -2 1 <(cat $out_dir/VIRSorter_global-phage-signal.csv | grep -v "^#" | cut -d ',' -f1| sort )  <(cat $out_dir/Metric_files/VIRSorter_affi-contigs.tab | grep -v '>' | tr '|' '\t' | awk '{print $1,"\t",$6}' | awk '{if($2!="-" && $2!="")print $0}' | sed -e 's/\-gene_[0-9]*//g' | sort -k 1,1 -k2,2 | uniq) \
		| tr ' ' '\t' | sort -k2 | uniq | cut -f2 | uniq -c | sort -r -k 1 -n | awk '{if($1!="1")print $2}' | sort >> /dev/shm/shared_phage_genes.txt

	##now we only want the phages that have at least one gene in the shared_phage_genes
	join -1 1 -2 1 <(cat $out_dir/VIRSorter_global-phage-signal.csv | grep -v "^#" | cut -d ',' -f1| sort )  <(cat $out_dir/Metric_files/VIRSorter_affi-contigs.tab | grep -v '>' | tr '|' '\t' | awk '{print $1,"\t",$6}' | awk '{if($2!="-" && $2!="")print $0}' | sed -e 's/\-gene_[0-9]*//g' | sort -k1,1 -k2,2 | uniq) \
		>> /dev/shm/phages_plus_phage_genes.txt

	sort -k2 /dev/shm/phages_plus_phage_genes.txt >> /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt

	##phages with shared gene content (pangenome)
	join -1 1 -2 2 -o2.1,1.1 /dev/shm/shared_phage_genes.txt /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt \
		| tr " " "\t" | cut -f1| sort | uniq >> /dev/shm/phages_with_shared_pangenome.txt

	echo "---------------------------------------------------"
	echo "Number of phages that have some shared gene content"
	cat /dev/shm/phages_with_shared_pangenome.txt | wc -l
}


## do not call this function on large concatenations, it is computationally expensive and very inefficient
function build_gene_presence_absence_matrix {

	out_dir=$1

	phages=`cat /dev/shm/phages_with_shared_pangenome.txt`
	genes=`cat /dev/shm/shared_phage_genes.txt`

	echo "starting building gene presence/absence matrix"
	##let us just take the unique sequences and the unique genes and then
	##loop over the sequences
	##    write the sequence name
	##    loop over the genes
	##        if found (grep if there is a line in the original file = very inefficient!)
	##            write 1
	##        else
	##            write 0
	rm -f /dev/shm/gene_presence_absence.txt

	line="sequence";for gene in $genes; do line="$line;$gene"; done ; echo $line >> /dev/shm/gene_presence_absence.txt

	##this can be computionally expensive because this is a [phage X gene] problem, and it greps with a regex. eg 1000 phages and 1000 genes will do 10^6 regex greps.
	line=""; for phage in $phages; do line="$phage" ; phage_esc=`echo $phage | sed -r 's/\[/\\\[/g' | sed -r 's/\]/\\\]/g'`; for gene in $genes; do gene_esc=`echo $gene | sed -r 's/\[/\\\[/g' | sed -r 's/\]/\\\]/g'`; count=`grep -c "$phage_esc $gene_esc" /dev/shm/phages_plus_phage_genes.txt`; line="$line;$count"; done; echo $line; done \
		>> /dev/shm/gene_presence_absence.txt

	echo "-----------------------------------------------------"
	echo "Gene presence/absence is in gene_presence_absence.txt"

	##copy to extra output file
	rm -f gene_presence_absence.txt
	cp /dev/shm/gene_presence_absence.txt gene_presence_absence.txt

}

function view_all_phages {
	out_dir=$1
	cat $out_dir/VIRSorter_global-phage-signal.csv | grep -v "^#" | cut -d ',' -f1| sort
}

function view_affi_all_contigs {
	##TODO: There is a bug here in the replacement of genes (it should only do this on the first column, not on the gene with the contig-name+gene
	##TODO: Improve performance (takes about 5 seconds on 1000 genomes)
	out_dir=$1

	cat $out_dir/Metric_files/VIRSorter_affi-contigs.tab | grep -v '>' | tr '|' '\t' \
	| awk '{if($6!="-" && $6!="")print $1,"\t",$6}' \
	| sed -e 's/\-gene_[0-9]*//g' | sed -e's/\./_/g' | sort -k 1,1 -k2,2 | uniq
}

function make_dictionary_of_annotations {

	rm -f /dev/shm/gi_dictionary
	rm -f /dev/shm/pc_dictionary
	rm -f /dev/shm/annotation_dictionary

	##individual gene information

	grep ">gi" virsorter-data/Phage_gene_catalog/phage_protein_14-03_RefseqABVir.faa | sed -e 's/>//g' | tr "." "_" \
		| awk -F "|" '{print $1"_"$2"_"$3"_"$4"_|"$5}' | sort | uniq >> /dev/shm/gi_dictionary


	grep "Phage_c" virsorter-data/Phage_gene_catalog/Phage_Clusters_current.tab | awk -F '|' '{print $1"|"$3"|"$4}' \
		| sed -e 's/\(Phage_cluster_\)\([0-9]*\)/\1\2\|\2/g' | sort -t '|' -k2,2 -n | cut -d '|' -f1,3,4 >> /dev/shm/pc_dictionary

	cat /dev/shm/gi_dictionary /dev/shm/pc_dictionary >> /dev/shm/annotation_dictionary
}

function view_top_annotations {
	##TODO: add extra fields, join does not work well now because of all the spaces


	rm -f /dev/shm/counts_gene_occurence_per_phage
	## join annotation directory
	## /dev/shm/annotation_dictionary

        cat /dev/shm/gene_occurence_per_phage | tr ' ' '\t' | sort -k2 | cut -f2 | uniq -c | sort -r -k 1 -n | sed -e 's/^\s*//g' | tr "." "_" | tr " " "|" \
		| sort -t '|' -k2 >> /dev/shm/counts_gene_occurence_per_phage

	##sort it back on occurence _after_ the join
	##you can also add -o2.3 if you like
	join -1 2 -2 1 -t "|" -o1.1,0,2.2 -a1 /dev/shm/counts_gene_occurence_per_phage /dev/shm/annotation_dictionary | uniq | sort -r -k 1 -n
}


function show_virsorter_summary {

	echo "summary of " $1
	out_dir=$1
	nr_top_rankings=10

	echo "Total number of viral sequences:"
	grep -c VIRSorter $out_dir/VIRSorter_global-phage-signal.csv

	echo "Number of entirely viral sequences:"
	grep ">" $out_dir/Predicted_viral_sequences/VIRSorter_*cat-*.fasta | \
	grep -E "cat_[1-3]" | wc -l

	echo "Number of prophages: "
	grep ">" $out_dir/Predicted_viral_sequences/VIRSorter_*cat-*.fasta | grep -E "cat_[4-6]" | wc -l

	show_iteration_summaries $out_dir

	echo "Final number of individual genes found:"
	cat $out_dir/Metric_files/VIRSorter_affi-contigs.tab | cut -d "|" -f6 | grep . | grep -v '-' | grep -c gi

	echo "Final number of gene clusters found:"
	cat $out_dir/Metric_files/VIRSorter_affi-contigs.tab | cut -d "|" -f6 | grep . | grep -v '-' | grep -c Phage_cl

	echo "-----------------------------------------------"
	echo "Final top individual genes and gene clusters found (all, not just the shared ones)"
	cat $out_dir/Metric_files/VIRSorter_affi-contigs.tab | cut -d "|" -f6 | grep . | grep -v '-' | sort | uniq -c | sort -r -k 1 -n | head -$nr_top_rankings



	##Only determine gene occurence per phage once (it contains a sort that takes multiple seconds) and store in memory:
	rm -f /dev/shm/gene_occurence_per_phage
	join -1 1 -2 1 <(view_all_phages $out_dir)  <(view_affi_all_contigs $out_dir) >> /dev/shm/gene_occurence_per_phage

	echo "-----------------------------------------------"
	echo "Top abundance of gene(cluster)s in pangenome"
	cat /dev/shm/gene_occurence_per_phage | tr ' ' '\t' | sort -k2 | uniq | cut -f2 | uniq -c | sort -r -k 1 -n \
		| awk '{if($1!="1")print $0}' | head -$nr_top_rankings

	echo "------------------------------------------------------------------"
	echo "Number of genes that occur at least in two of the predicted phages"
	cat /dev/shm/gene_occurence_per_phage | tr ' ' '\t' | sort -k2 | uniq | cut -f2 | uniq -c | sort -r -k 1 -n | awk '{if($1!="1")print $0}' | wc -l


	##determine pangenome
	determine_shared_genes $out_dir

	##WATCH OUT BEFORE OUTCOMMENTING THIS ONE: COMPUTATIONALLY VERY EXPENSIVE!!
	## build_gene_presence_absence_matrix $out_dir

	echo "---------------------------------------------------"
	echo "VirSorter annotation of top 20 phage clusters found"
	echo "---------------------------------------------------"

	make_dictionary_of_annotations

	view_top_annotations | head -20

	##TODO 
	##	(1) prepare a dictionary for both the individual genes and the phage clusters so that you can join this with the top occurences of both of them
	##	    see above "Top abundance of gene (cluster)s in pangenome"
	##	(2) prevent code repetition by storing subresults of data sets or put these in functions
	##	(3) do not grep in a while loop but also use a join for the last combination with the annotation information (with the prepared dictionary)
	##

	## the next line is also a possibility but there is no outer join
	##cat /dev/shm/counts_gene_occurence_per_phage | sort -r -k 1 -n | grep -v "^1|" | awk -F '|' '{print $2}' \
	##	| while read LINE; do grep "$LINE|" /dev/shm/annotation_dictionary; done | cut -d '|' -f1,2,3


	echo "--------------------------------------------------------------------"
	echo "VirSorter annotation of top 20 phage clusters found within pangenome"
	echo "--------------------------------------------------------------------"
	cat /dev/shm/gene_occurence_per_phage \
		| tr ' ' '\t' | sort -k2 | uniq | cut -f2 | uniq -c | sort -r -k 1 -n | awk '{if($1!="1")print $0}' | awk '{print $2}' \
		| grep Phage_cluster | while read LINE; do grep "$LINE|" virsorter-data/Phage_gene_catalog/Phage_Clusters_current.tab; done | tr '|' '\t' | cut -f1,3 | head -20

	echo "-----------------------------------------------------------------------"
	echo "VirSorter annotation of top 20 individual phage genes within pangenome" 
	echo "-----------------------------------------------------------------------"

	cat /dev/shm/gene_occurence_per_phage \
		| tr ' ' '\t' | sort -k2 | uniq | cut -f2 | uniq -c | sort -r -k 1 -n | awk '{if($1!="1")print $0}' | awk '{print $2}' \
		| grep gi_ | tr "." "_" | while read gene; do grep $gene /dev/shm/gi_dictionary ; done | head -20

	echo "------------------------------------------"
	echo "Top most shared genes and gene clusters"
	cat $out_dir/Metric_files/VIRSorter_affi-contigs.tab | grep -v '>' | tr '|' '\t' | awk '{print $1"\t"$6}' | awk '{if($2!="-" && $2!="")print $0}' | sed -e 's/\-gene_[0-9]*//g' | sort -k 2,2 -k1,1 | cut -f2 | uniq -c | sort -r -k 1 -n | head -$nr_top_rankings
}

## this is just an alias (do not add any programming logic)
function show_vs {
	out_dir=$1
	show_virsorter_summary $out_dir
}
