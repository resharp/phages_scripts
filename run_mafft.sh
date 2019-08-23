
#example workflow
#	run_mafft PC_1				#multiple alignment, store in *_mafft.fasta file
#	build_profile PC_1			#build hmm profile, store in *.hmm file
#	search_own_cluster_with_profile PC_1	#search with profile, store in *hmm_results.txt file and table file
#	search_all_proteins_with_profile PC_1	#search with profile, for all proteins!

function run_mafft {

	PC=$1

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	fasta_file=$gene_dir/mcl_75.$sample/${PC}.fasta
	out_file=$gene_dir/mcl_75.$sample/${PC}_mafft.fasta

	echo $fasta_file
	echo $out_file

	mafft --quiet $fasta_file > $out_file
}

function cleanup_mafft {

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	out_files=$gene_dir/mcl_75.$sample/*mafft.fasta

	#remove all out files
	rm -f $out_files
}


#TODO: this function calculate_gaps is not the way to go
#use trimall to trim the alignment, and then determine what you want alignments are still useful
function calculate_gaps {

	PC=$1

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	out_file=$gene_dir/mcl_75.$sample/${PC}_mafft.fasta

	echo "Nr of gaps and alignments protein cluster: " $PC

	nr_gaps=$(cat $out_file | grep -v ">" | tr -d '\n' | grep -o -i '-' | wc -l)
	nr_alignments=$(cat $out_file | grep -v ">" | tr -d '\n' | grep -o -i '[A-Z]' | wc -l)

	nr_total=$((nr_gaps+nr_alignments))

	echo $nr_gaps $nr_total | awk '{print $1/$2}'

}

function build_profile {

	PC=$1

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	align_file=$gene_dir/mcl_75.$sample/${PC}_mafft.fasta
	hmm_file=$gene_dir/mcl_75.$sample/${PC}_mafft.hmm

	hmmbuild --cpu 12 $hmm_file $align_file

}


#this only searches the sequences in the cluster itself
function search_own_cluster_with_profile {

	PC=$1

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	fasta_file=$gene_dir/mcl_75.$sample/${PC}.fasta
	hmm_file=$gene_dir/mcl_75.$sample/${PC}_mafft.hmm
	out_file=$gene_dir/mcl_75.$sample/${PC}_mafft_hmm_results_cluster.txt
	out_table=$gene_dir/mcl_75.$sample/${PC}_mafft_hmm_results_cluster_table.txt

	echo $fasta_file
	echo $hmm_file

	#search all proteins in cluster with the hmm profile of the cluster
	hmmsearch --cpu 12 --tblout $out_table $hmm_file $fasta_file > $out_file 

}

function search_all_proteins_with_profile {

	PC=$1

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	#NB fasta file with all samples!
	# * *  *   *    *
	fasta_file=$gene_dir/gene_samples_simple.fasta
	hmm_file=$gene_dir/mcl_75.$sample/${PC}_mafft.hmm
	out_file=$gene_dir/mcl_75.$sample/${PC}_mafft_hmm_results_all.txt
	out_table=$gene_dir/mcl_75.$sample/${PC}_mafft_hmm_results_all_table.txt

	echo $fasta_file
	echo $hmm_file
	echo $out_file
	echo "table output in: " $out_table

	#search all proteins with hmm profile of one cluster
	hmmsearch --tblout $out_table $hmm_file $fasta_file > $out_file

}


function example_compare_hmm_hits {

	paste 	<(cat /hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000/mcl_75.I20/PC_1_mafft_hmm_results_table.txt | grep -v "#" | awk '{print $1}')\
		<(cat /hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000/mcl_75.I20/PC_1_mafft_hmm_results_all_table.txt | grep -v "#" | awk '{print $1}') | less
}


function run_alignment_evaluation_workflow {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	sample="I20"
	nr_pcs=$(wc -l $gene_dir/out.pw_blastout_mcl_75.abc.$sample | awk '{print $1}')

	log_file=$gene_dir/log_evaluation.pw_blastout_mcl_75.abc.$sample.txt

	echo "alignments, hmm profiles and profile searches for all clusters in "
	echo "total nr of protein clusters" $nr_pcs

	#NB: this operates on all pcs, when testing set head!
	# * *  *   *     *
	#pcs=$(cat --number $gene_dir/out.pw_blastout_mcl_75.abc.$sample | awk '{print "PC_" $1}' | head -10 )
	pcs=$(cat --number $gene_dir/out.pw_blastout_mcl_75.abc.$sample | awk '{print "PC_" $1}')

	#empty log_file
	> $log_file
	for pc in $pcs
	do
		echo "running workflow for cluster "$pc
		run_mafft $pc >> $log_file
		build_profile $pc >> $log_file
		search_all_proteins_with_profile $pc >> $log_file
	done
}
