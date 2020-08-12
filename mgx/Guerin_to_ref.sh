
# first cd to /hosts/linuxhome/chaperone/richard/guerin_data
# sample_dir=/hosts/linuxhome/chaperone/richard/guerin_data
# ref_file=~/scripts/mgx/ref_seqs/ids_ref_genomes.txt

# cd $sample_dir


# Guerin.fna downloaded from https://github.com/linsalrob/crAssphage/blob/master/Guerin_Phages/Guerin.fna
# 
# split Guerin.fna in different fasta files in different subdirectories for all reference genomes
#
function split_fasta_by_ids {

	ref_file=$1

        if [ -z "$ref_file" ]
        then
                echo "please provide ref_file as argument (file with ref genome ids)"
        else

		conda deactivate
		conda activate python37
		#seqtk is in env python37
		seqtk subseq Guerin.fna ${ref_file} > Guerin_subset.fa

		samples=$(grep ">" Guerin_subset.fa | cut -d " " -f1 | sed -e 's/>//g')
		# for sample in $samples; do echo $sample; done
		for sample in $samples
			do
				rm -f -r ${sample}
				mkdir ${sample}
				grep -A1 -w $sample Guerin_subset.fa > ${sample}/${sample}.fasta
			done
	fi
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

		echo $sample

		# two sequences (for candidate genus 7 and 8) use a different genetic code: genetic code 15
		# use prodigal -g 15
		# inf125_s_2
		# eld241-t0_s_1
		if [[ $sample == "inf125_s_2" || $sample == "eld241-t0_s_1" ]]; then
			echo "use different genetic code 15"
			# NB: Don't use meta mode if the translation table is provided! see:
			# https://github.com/merenlab/anvio/issues/1074
			prodigal -i $file -g 15 -o ${sample}/${sample}.genes -a ${sample}/${sample}.proteins.faa
		else
			echo "use default genetic code 11"
			prodigal -i $file -o ${sample}/${sample}.genes -a ${sample}/${sample}.proteins.faa -p meta
		fi

	done

}
