#
# run_crassphage_samples
#
#---------------
source_dir=/hosts/linuxhome/mgx/DB/MGXDB/MGXDB

#TODO: replace and set at location on mutant
sample_dir=/hosts/linuxhome/chaperone/tmp/richard/crassphage_samples

master_dir=/hosts/linuxhome/chaperone/tmp/richard/crassphage_samples

function run_crassphage_samples {

	copy_reference_genome

	select_samples

	run_alignment

}

function copy_reference_genome {

	#used by select_samples: directories of 18,417 samples
	cp ${master_dir}/MGXDB_samples.txt ${sample_dir}

	#reference genome + index for bwa mem
	cp ${master_dir}/crassphage* ${sample_dir}

}


#copy some samples from
# /hosts/linuxhome/mgx/DB/MGXDB/MGXDB
function select_samples {

	samples=$(cat $sample_dir/MGXDB_samples.txt | head -200 | awk '{print $9}')

	nr_selected=0
	for sample in $samples
	do
		# echo $source_dir/$sample/${sample}_metadata
		#todo: does the metadata file exist?
		nr_gut=0
		nr_dorei=0
		meta_file=$source_dir/$sample/${sample}_metadata
		if [[ -f "$meta_file" ]]; then
			nr_gut=$(grep -c "gut" $source_dir/$sample/${sample}_metadata)
			nr_dorei=$(grep -c "dorei" $source_dir/$sample/${sample}_metadata)
		fi
		if [[ $nr_dorei == "1" && $nr_gut == "1" ]]; then
			echo "'dorei' and 'gut' found in" $sample
			gz_file=$source_dir/$sample/${sample}_filtered.fastq.gz
			if [ -f "$gz_file" ];then
				echo $gz_file
				nr_selected=$((nr_selected+1))
				cp $gz_file $sample_dir/${sample}_filtered.fastq.gz
				echo $sample >> $sample_dir/MGXDB_selected_samples.txt
			fi
		fi
	done
	echo $nr_selected
}


#what is the size of the reads?
#optional: downsample the samples, e.g. on 10,000 short reads

#run alignment (try out multiple aligners)
#run bwa
function run_alignment {

	conda activate mgx

	samples=$(cat $sample_dir/MGXDB_selected_samples.txt)

	for sample in $samples; do

		file=$sample_dir/${sample}_filtered

		gunzip ${file}.fastq.gz

		bwa mem $sample_dir/crassphage_refseq.fasta ${file}.fastq > ${file}.sam

		samtools view -S -b $file.sam > $file.bam

		samtools sort -o $file.sorted.bam -@ 4 $file.bam

		samtools index $file.sorted.bam

		samtools idxstats $file.sorted.bam > $file.sorted.idstats.txt
	done
}

#conda activate samtools
#inspect stats
