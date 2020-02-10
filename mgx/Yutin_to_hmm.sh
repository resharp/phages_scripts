
gene_dir=/hosts/linuxhome/chaperone/tmp/richard/yutin_data/msa
ref_dir=/hosts/linuxhome/chaperone/tmp/richard/guerin_data


function download_alignments {

	#We downloaded the multiple sequence alignments for building HMM profiles for every annotated protein family.
	#See directory msa
	#(files in there have a afa extension "aligned fasta")

	#See README on
	#ftp://ftp.ncbi.nih.gov/pub/yutinn/crassphage_2017/

	#example only get protein alignments:
	wget --recursive --no-parent --cut-dirs=4 ftp://ftp.ncbi.nih.gov/pub/yutinn/crassphage_2017/prot_align/
}



function build_profiles {

	conda deactivate
	conda activate mgx

	gene_dir=$1

	files=$(find $gene_dir/*.afa)
	for align_file in $files
	do
		# echo $file
		base_file=$(basename "$align_file")

		gene=${base_file%.afa}
		echo $gene

		hmm_file=${gene_dir}/${gene}.hmm

		echo $align_file
		echo $hmm_file
		hmmbuild --cpu 12 $hmm_file $align_file
	done
}


function search_genomes {


	#now we want to search all ref genomes for all 
	ref_dir=$1
	gene_dir=$2

	profiles=$(find $gene_dir/*.hmm)

	subdirs=$(find $ref_dir -type d ! -name "guerin_data")

	for subdir in $subdirs
	do
		genome=$(basename "$subdir")

		fasta_file=${ref_dir}/${genome}/${genome}.proteins.faa

		for profile in $profiles
		do
			echo "one combination"

			base_file=$(basename "$profile")
			gene=${base_file%.hmm}
			echo $gene

			out_file=${ref_dir}/${genome}/${genome}-${gene}.txt
			out_table=${ref_dir}/${genome}/${genome}-${gene}_table.txt

			#echo $out_file
			#echo $out_table

			# search all genes in every genome
			hmmsearch --cpu 12 --tblout $out_table $profile $fasta_file > $out_file
		done
	done
}


# show all hmm hits on all ref genomes
function show_profile_hits {

	ref_dir=/hosts/linuxhome/chaperone/tmp/richard/guerin_data

	# last column AA length calculation: diff between end and start pos of CDS, substract 3 (stop codon) and divide by 3 for AA length
	cat $ref_dir/*/*_table.txt | grep -v "^#" | awk '{print $1"\t"$3"\t"$5"\t"$6"\t"$20"\t"$22"\t"($22-$20+1-3)/3}' | sort -k2,2 -k1,1 | uniq | less -x 30

}
