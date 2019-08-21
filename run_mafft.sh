#running mafft
function run_mafft {

	PC=$1

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	fasta_file=$gene_dir/mcl.$sample/${PC}.fasta
	out_file=$gene_dir/mcl.$sample/${PC}_mafft.fasta

	echo $fasta_file
	echo $out_file

	mafft $fasta_file > $out_file
}

function cleanup_mafft {

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	out_files=$gene_dir/mcl.$sample/*mafft.fasta

	#remove all out files
	rm -f $out_files
}


#TODO: this function calculate_gaps is not the way to go
#use trimall to trim the alignment, and then determine what you want alignments are still useful
function calculate_gaps {

	PC=$1

	sample="I20"

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	out_file=$gene_dir/mcl.$sample/${PC}_mafft.fasta

	echo "Nr of gaps and alignments protein cluster: " $PC

	nr_gaps=$(cat $out_file | grep -v ">" | tr -d '\n' | grep -o -i '-' | wc -l)
	nr_alignments=$(cat $out_file | grep -v ">" | tr -d '\n' | grep -o -i '[A-Z]' | wc -l)

	nr_total=$((nr_gaps+nr_alignments))

	echo $nr_gaps $nr_total | awk '{print $1/$2}'

}
