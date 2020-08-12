# run pipeline on one SRA project

# e.g. run_all_samples $source_dir $sample_dir 10

# make list of samples

# run complete pipeline on those samples
# copy to sample_dir

function run_all_samples {

	source_dir=$1
	sample_dir=$2
	nr_samples=$3

	files_1=$(find $source_dir/*_1.fastq.gz | head -${nr_samples})

        for file_1 in $files_1
        do
                base_file_1=$(basename "$file_1")
                # e.g. run=ERR525689

                run=${base_file_1%_1.fastq.gz}

		run_sample $source_dir $sample_dir $run
        done

}

function run_sample {

        source_dir=$1
        sample_dir=$2
        run=$3

	copy_and_run_trimming $source_dir $sample_dir $run

	run_mapping $sample_dir $run

	nr_mapped_reads=$(head -1 ${sample_dir}/${run}/${run}.sorted.idstats.txt | cut -f3)

	if [ $nr_mapped_reads -gt 0 ]
	then
		echo "Nr of reads mapped: "$nr_mapped_reads
		run_diversiutils $sample_dir $run
		run_calc_measures $sample_dir $run
	else
		echo "No reads mapped to reference genome!"
	fi
}


# debug_sample $source_dir $sample_dir ERR1136746
function debug_sample {

	conda deactivate
	conda activate mgx

        source_dir=$1
        sample_dir=$2
        run=$3

	copy_no_trimming $source_dir $sample_dir $run

	run_mapping $sample_dir $run

	nr_mapped_reads=$(head -1 ${sample_dir}/${run}/${run}.sorted.idstats.txt | cut -f3)

	if [ $nr_mapped_reads -gt 0 ]
	then
		echo "Nr of reads mapped: "$nr_mapped_reads
		run_diversiutils $sample_dir $run
		run_calc_measures $sample_dir $run
	else
		echo "No reads mapped to reference genome!"
	fi

}

function copy_no_trimming {
	source_dir=$1
	sample_dir=$2
	run=$3

	echo $source_dir
	echo $sample_dir
	echo $run

        file_1=${source_dir}/${run}_1.fastq.gz
        file_2=${source_dir}/${run}_2.fastq.gz

	mkdir $sample_dir/${run}
	# cp $file_1 $sample_dir/${basename}_1.fastq.gz
	/bin/cp -rfv $file_1 $sample_dir/${run}/
	/bin/cp -rfv $file_2 $sample_dir/${run}/
}

function run_mapping {

	sample_dir=$1
	run=$2
	file_1=${sample_dir}/${run}/${run}_1.fastq.gz
	file_2=${sample_dir}/${run}/${run}_2.fastq.gz

	ll $file_1
	ll $file_2

	#base for all other bam/same/etc files
	file=$sample_dir/${run}/${run}

	ref_seq=scripts/mgx/ref_seqs/crassphage_refseq.fasta

        #TODO: check bwa mem options
	# bwa mem should work on .gz files
        bwa mem -t 16 ${ref_seq} ${file_1} ${file_2} > ${file}.sam


	ll ${file}.sam
	#TODO: check if you can pipe output directly to prevent space usage
        #convert sam to bam
        samtools view -@ 16 -S -b $file.sam > $file.bam

	ll ${file}.bam
	#rm -f $file.sam


        #sort bam
        samtools sort -@ 16 -o $file.sorted.bam $file.bam

	ll ${file}.sorted.bam
	#rm -f $file.bam

        #index bam
        samtools index $file.sorted.bam

        #create stats for mapped reads
        samtools idxstats $file.sorted.bam > $file.sorted.idstats.txt
}

function run_diversiutils {

        sample_dir=$1
        sample=$2

	ref_seq=scripts/mgx/ref_seqs/crassphage_refseq.fasta

        tools/DiversiTools/bin/diversiutils_linux -bam $sample_dir/${sample}/${sample}.sorted.bam\
                -ref ${ref_seq}\
                -orfs scripts/mgx/crassphage_codingregions.txt\
                -stub $sample_dir/${sample}/${sample}

        #remove the "<NA>" strings
        sed -e 's/<NA>//g' $sample_dir/${sample}/${sample}_AA.txt > $sample_dir/${sample}/${sample}_AA_clean.txt

	rm -f $sample_dir/${sample}/${sample}_AA.txt
}

function run_calc_measures {

        sample_dir=$1
        sample=$2

        conda deactivate
        conda activate python37

	/bin/cp -rfv source/phages/codon_syn_non_syn_probabilities.txt ${sample_dir}/

        #create measures (based on ${sample}_AA_clean.txt)
        python source/phages/CalcDiversiMeasures.py -d $sample_dir -s $sample
}
