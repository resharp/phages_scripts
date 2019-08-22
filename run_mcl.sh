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

function run_mcl {

	#default inflation is 2.0. we run with 2.5 for comparing with cytoscape plug-in
	#keep it simple run with all default options

	#TODO: change location, this one is only running on 1000 genomes
	#---- * *  *     *               *
	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	#run with different inflations

	samples="1.5 2.0 2.5 3.0"
	for sample in $samples:
	do
		mcl $gene_dir/pw_blastout_mcl.abc --abc -I ${sample} -te 4 --d
	done

	samples="I15 I20 I25 I30"

	for sample in $samples
		do echo "largest clusters for inflation: "$sample
		cat $gene_dir/out.pw_blastout_mcl.abc.$sample | awk '{print NF" proteins"}' | head -3
	done

}

function split_mcl {

	#samples denote runs for different mcl inflation factors
	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	mkdir $gene_dir/mcl.$sample

	n=1
	cat $gene_dir/out.pw_blastout_mcl.abc.$sample |\
		while read line
		do
			for word in $line
			do
				echo $word >> $gene_dir/mcl.$sample/PC_$n.txt
			done
			n=$((n+1))
		done
}

function split_fasta_according_to_mcl {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	sample="I20"

	files=$(find $gene_dir/mcl.$sample/PC_*.txt)
	for file in $files
	do
		#run your tool here
		#https://github.com/lh3/seqtk

		fasta_file=$(echo $file | sed -e 's/\.txt/\.fasta/g' )

		echo "creating" $fasta_file
		seqtk subseq $gene_dir/gene_samples_simple.fasta $file > $fasta_file
	done
}
