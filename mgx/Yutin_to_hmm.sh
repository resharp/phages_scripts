
#gene_dir=/hosts/linuxhome/chaperone/tmp/richard/yutin_data/msa
#ref_dir=/hosts/linuxhome/chaperone/tmp/richard/guerin_data


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

	conda deactivate
	conda activate mgx

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

	# this result has been written to hmm_hits.tsv
}

function show_genes {

	ref_dir=/hosts/linuxhome/chaperone/tmp/richard/guerin_data

	# NB: cut -f8 gives the subdirectory name with the ref genome name
	paste	<(grep ">" ${ref_dir}/*/*.proteins.faa | cut -f8 -d "/") \
		<(grep ">" ${ref_dir}/*/*.proteins.faa | cut -f2 -d ">" | cut -f1 -d " ") \
		<(grep ">" ${ref_dir}/*/*.proteins.faa | cut -d ">" -f2 | cut -d "#" -f2-)

	# this result has been written to ref_genes.tsv
}

# also run a blastp search of all proteins against all proteins
# it may make it possible to link the proteins to annotated proteins of the p-crAssphage, for proteins not hit by a HMM profile
function blastp_all_against_all {

	ref_dir=/hosts/linuxhome/chaperone/tmp/richard/guerin_data

	# use this for building the
	cat ${ref_dir}/*/*.proteins.faa | awk '{ if($0~">") {print $1} else print $0  }' > ${ref_dir}/all_refs.proteins.faa

	makeblastdb -in ${ref_dir}/all_refs.proteins.faa -parse_seqids -title "ref genomes genes" -dbtype prot -out ${ref_dir}/all_refs.proteins.db.faa

	evalue=1e-3

	echo 'start running blastp'
	blastp -query $ref_dir/all_refs.proteins.faa -db $ref_dir/all_refs.proteins.db.faa\
		-evalue $evalue\
		-out $ref_dir/pairwise_blastout.txt\
		-num_threads 12\
		 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
	#-outfmt '7 qseqid sseqid length qlen qstart qend sstart send evalue'

	echo 'filtering out protein comparisons against themselves'
	cat $ref_dir/pairwise_blastout.txt | awk '{if ($1!=$2) print $0}' > $ref_dir/pairwise_blastout_filtered.txt

	echo 'number of hits:'
	wc -l $ref_dir/pairwise_blastout_filtered.txt

}
