
# first cd to
# /hosts/linuxhome/chaperone/tmp/richard/guerin_data
# sample_dir=/hosts/linuxhome/chaperone/tmp/richard/guerin_data
# cd $sample_dir

function split_fasta_by_ids {

	conda deactivate
	conda activate python37
	#seqtk is in env python37
	seqtk subseq Guerin.fna ids_ref_genomes.txt > Guerin_subset.fa

	samples=$(grep ">" Guerin_subset.fa | cut -d " " -f1 | sed -e 's/>//g')
	# for sample in $samples; do echo $sample; done
	for sample in $samples
		do
			rm -f -r ${sample}
			mkdir ${sample}
			grep -A1 -w $sample Guerin_subset.fa > ${sample}/${sample}.fasta
		done
}


function run_prodigal_on_ref {


	conda deactivate
	conda activate prodigal

	# prodigal -i $fasta_file -o annotations/my.genes -a annotations/my.proteins.faa -p meta
	files=$(find */*.fasta)
	for file in $files
	do
		name=$(basename $file)
		sample=$(echo $name | cut -f1 -d '.')

		# echo $sample
		# here we should run it on a single genomen (does it matter?)
		prodigal -i $file -o ${sample}/${sample}.genes -a ${sample}/${sample}.proteins.faa -p meta

	done
	# prodigal -i $file -o annotations/my.genes -a annotations/my.proteins.faa -p meta
}
