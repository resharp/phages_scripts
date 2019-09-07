#example mcl workflow with preprocessing and postprocessing steps
# depends on functions in view_stats.sh!
function run_mcl_workflow {

	#TODO: change location, this one is running on ALL genomes
	#---- * *  *     *               *
	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_all
	#gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	#depends on result op pairwise blastp
	prepare_abc_for_mcl75 $gene_dir

	conda deactivate
	conda activate virsorter
	#note: change the required inflation factors  in script
	run_mcl $gene_dir

	#build files for evaluating mcl clustering with hmm profiles (depends on seqtk in conda env python37)
	#resulting files are input for MclClusterEvaluation.py
	conda deactivate
	conda activate python37

	#NB: add extra lines for different evaluations I15, I20, I25, etc ..
	split_mcl $gene_dir I25
	split_fasta_according_to_mcl $gene_dir I25

	#build files for calculating shared gene content (dependency on functions in view_stats.sh)
	#resulting files are input for PhageSharedContent.py
	make_phage_ip_table $gene_dir
	make_phage_ip_table_short $gene_dir
	make_phage_table $gene_dir

	#also builds pc_table
	make_ip_pc_table $gene_dir I25
	make_filtered_ip_pc_table $gene_dir "mcl_75.I25" 100

}

function prepare_sif_for_cytoscape {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes

	#http://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html
	#sif format example:
	#node1 typeA node2
	#node2 typeB node3 node4 node5
	cat $gene_dir/pairwise_blastout_filtered.txt | awk '{print $1 " homolog " $2}' > $gene_dir/pw_blastout_cytoscape.sif
}

function prepare_abc_for_mcl {

	gene_dir=$1

	#prepare format for mcl (without weight)
	cat $gene_dir/pairwise_blastout_filtered.txt | awk '{print $1 "\t" $2}' > $gene_dir/pw_blastout_mcl.abc

}

function prepare_abc_for_mcl75 {


	gene_dir=$1

	#prepare format for mcl (without weight)
	cat $gene_dir/pairwise_blastout_filtered_75coverage.txt | awk '{print $1 "\t" $2}' > $gene_dir/pw_blastout_mcl_75.abc

}

#mcl is installed in conda env virsorter
function run_mcl {

	#default inflation is 2.0. we run with 2.5 for comparing with cytoscape plug-in
	#keep it simple run with all default options

	#TODO: change location, this one is running on all genomes
	#---- * *  *     *               *
	gene_dir=$1

	file_in_mcl=pw_blastout_mcl_75.abc
	#run with different inflations

	#samples="1.5 2.0 2.5 3.0"
	samples="2.5"
	for sample in $samples:
	do
		mcl $gene_dir/$file_in_mcl --abc -I ${sample} -te 4 --d
	done

	#samples="I15 I20 I25 I30"
	samples="I25"

	for sample in $samples
		do echo "largest clusters for inflation: "$sample
		cat $gene_dir/out.$file_in_mcl.$sample | awk '{print NF" proteins"}' | head -3
	done

}

#example split_mcl $gene_dir I25
# builds set of files in directory (file for each PC consisting of list of IPs)
# this is used by split_fasta_according_to_mcl
#	mcl_75.$sample/PC_$n.txt
function split_mcl {

	gene_dir=$1

	#samples denote runs for different mcl inflation factors
	sample=$2

	##NB processes 75 coverage file
	# * *  *   *    *
	file_out_mcl=out.pw_blastout_mcl_75.abc

	#mkdir $gene_dir/mcl.$sample
	mkdir $gene_dir/mcl_75.$sample

	n=1
	cat $gene_dir/$file_out_mcl.$sample |\
		while read line
		do
			for word in $line
			do
				echo $word >> $gene_dir/mcl_75.$sample/PC_$n.txt
			done
			n=$((n+1))
		done
}

#uses seqtk in python37 conda env
#example  split_fasta_according_to_mcl $gene_dir I25
#
# depends on 
#	pc_table.mcl_75.I25 or .$sample
#	the PC_1.txt files etc in the  mcl_75.I25 directory, or mcl_75.$sample (built with split_mcl)
function split_fasta_according_to_mcl {

	gene_dir=$1

	sample=$2

	pcs=$(cat $gene_dir/pc_table.mcl_75.$sample)

	for pc in $pcs
	do
		fasta_file=$gene_dir/mcl_75.$sample/$pc.fasta
		file=$gene_dir/mcl_75.$sample/$pc.txt

		echo "creating" $fasta_file
		seqtk subseq $gene_dir/gene_samples_simple.fasta $file > $fasta_file
	done
}
