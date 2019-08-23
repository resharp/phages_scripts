function prepare_sif_for_cytoscape {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes

	#http://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html
	#sif format example:
	#node1 typeA node2
	#node2 typeB node3 node4 node5
	cat $gene_dir/pairwise_blastout_filtered.txt | awk '{print $1 " homolog " $2}' > $gene_dir/pw_blastout_cytoscape.sif
}

function prepare_abc_for_mcl {


	#TODO: change location, this one is only running on 1000 genomes
	#---- * *  *     *               *
	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	#prepare format for mcl (without weight)
	cat $gene_dir/pairwise_blastout_filtered.txt | awk '{print $1 "\t" $2}' > $gene_dir/pw_blastout_mcl.abc

}

function prepare_abc_for_mcl75 {


	#TODO: change location, this one is only running on 1000 genomes
	#---- * *  *     *               *
	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	#prepare format for mcl (without weight)
	cat $gene_dir/pairwise_blastout_filtered_75coverage.txt | awk '{print $1 "\t" $2}' > $gene_dir/pw_blastout_mcl_75.abc

}

#mcl is installed in conda env virsorter
function run_mcl {

	#default inflation is 2.0. we run with 2.5 for comparing with cytoscape plug-in
	#keep it simple run with all default options

	#TODO: change location, this one is only running on 1000 genomes
	#---- * *  *     *               *
	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	file_in_mcl=pw_blastout_mcl_75.abc
	#run with different inflations

	samples="1.5 2.0 2.5 3.0"
	for sample in $samples:
	do
		mcl $gene_dir/$file_in_mcl --abc -I ${sample} -te 4 --d
	done

	samples="I15 I20 I25 I30"

	for sample in $samples
		do echo "largest clusters for inflation: "$sample
		cat $gene_dir/out.$file_in_mcl.$sample | awk '{print NF" proteins"}' | head -3
	done

}

#example split_mcl I25
function split_mcl {

	#samples denote runs for different mcl inflation factors
	sample=$1

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

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
#example  split_fasta_according_to_mcl I25
function split_fasta_according_to_mcl {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	sample=$1

	files=$(find $gene_dir/mcl_75.$sample/PC_*.txt)

	for file in $files
	do
		#run your tool here
		#https://github.com/lh3/seqtk

		fasta_file=$(echo $file | sed -e 's/\.txt/\.fasta/g' )

		echo "creating" $fasta_file
		seqtk subseq $gene_dir/gene_samples_simple.fasta $file > $fasta_file
	done
}
