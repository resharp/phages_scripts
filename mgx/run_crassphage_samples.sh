#
# run_crassphage_samples
#
#---------------
admin_dir=/hosts/linuxhome/mgx/DB/MGXDB
source_dir=/hosts/linuxhome/mgx/DB/MGXDB/MGXDB

#TODO: replace and set at right location on mutant
sample_dir=/hosts/linuxhome/mutant31/tmp/richard/crassphage_samples
#sample_dir=/hosts/linuxhome/chaperone/tmp/richard/crassphage_samples
master_dir=/hosts/linuxhome/chaperone/tmp/richard/crassphage_samples

function copy_reference_genome {

	#used by select_samples: directories of 18,417 samples
	cp ${master_dir}/MGXDB_samples.txt ${sample_dir}

	#reference genome + index for bwa mem
	cp ${master_dir}/crassphage* ${sample_dir}

}

#this is obsolete, see Python
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

#this is obsolete, see Python
function select_crassphage {

	samples=$(cat $sample_dir/MGXDB_samples.txt | awk '{print $9}')

	> $sample_dir/MGXDB_crassphage.txt
	nr_selected=0
	for sample in $samples
	do
		nr_crass=0
		meta_file=$source_dir/$sample/${sample}_metadata
		if [[ -f "$meta_file" ]]; then
			nr_crass=$(grep -c "crAssphage" $source_dir/$sample/${sample}_metadata)
		fi
		if [[ $nr_crass == "1" ]]; then
			echo "'crAsshage' found in" $sample
			gz_file=$source_dir/$sample/${sample}_filtered.fastq.gz
			if [ -f "$gz_file" ];then
				echo $gz_file
				nr_selected=$((nr_selected+1))
				#cp $gz_file $sample_dir/${sample}_filtered.fastq.gz
				echo $sample >> $sample_dir/MGXDB_crassphage.txt
			fi
		fi
	done
	echo "Nr containing assembled crassPhage: "$nr_selected
}



#what is the size of the reads?
#optional: downsample the samples, e.g. on 10,000 short reads

#run alignment (try out multiple aligners)
#run bwa
function run_crassphage_samples {

	conda deactivate
	conda activate mgx

	samples=$(cat $admin_dir/taxon_1211417_counts_sorted.txt | awk '{if ($3>10000) print $1}' | grep -v "sample")

	for sample in $samples; do
		echo $sample
		#run_sample $sample
	done
}

function run_one_sample {

	conda deactivate
	conda activate mgx

	sample=MGXDB000864 #this is the last of the > 10,000 total_count, so we have to exclude this in the final sample
	run_sample $sample
}

function run_sample {

	sample=$1

	echo "running sample " $sample

	#make a directory for each sample to copy the gz file to
	mkdir $sample_dir/${sample}

	#TODO
	#copy gz file from source_dir
	echo "copying gz file"
	gz_file=$source_dir/$sample/${sample}_filtered.fastq.gz
	if [ -f "$gz_file" ];then

		echo $gz_file

		rsync -av --force $gz_file $sample_dir/${sample}
	fi

	file=$sample_dir/${sample}/${sample}_filtered

	# to do: no need to unzip for bwa mem
	echo "unzip gz file"
	gunzip ${file}.fastq.gz

	#TODO: check bwa mem options
	bwa mem -t 16 $sample_dir/crassphage_refseq.fasta ${file}.fastq > ${file}.sam

	#convert sam to bam
	samtools view -@ 16 -S -b $file.sam > $file.bam

	#sort bam
	samtools sort -@ 16 -o $file.sorted.bam $file.bam

	#index bam
	samtools index $file.sorted.bam

	#create stats for mapped reads
	samtools idxstats $file.sorted.bam > $file.sorted.idstats.txt

	#remove intermediary files (only keep the sorted.bam file)
	rm -f ${file}.fastq
	rm -f ${file}.sam
	rm -f ${file}.bam
	#TODO: if we integrate the diversitools post processing step, we can also remove the sorted.bam file

	#copy to chaperone
	rsync -av --force $sample_dir/${sample} $master_dir

	#remove sorted.bam
	rm -f ${file}.sorted.bam
	rm -f ${file}.sorted.bam.bai

}

#conda activate samtools
#inspect stats
