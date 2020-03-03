# run complete pipeline on one SRA project
# for more than one reference genome
# e.g.

# e.g. run_all_samples $source_dir $sample_dir 10
#source_dir=/hosts/linuxhome/mutant14/tmp/richard/sra/ERP005989
#sample_dir=/hosts/linuxhome/mutant26/tmp/richard/ERP005989


# this contains one sample that definitily should have reads that are supposed to map against the CRassphage because it does
# when processing the files from our MGXDB database earlier
#source_dir=/hosts/linuxhome/mutant14/tmp/richard/sra/PRJEB11532
#sample_dir=/hosts/linuxhome/mutant26/tmp/richard/PRJEB11532

# run complete pipeline on those samples
function run_all_samples {

	source_dir=$1
	sample_dir=$2
	nr_samples=$3
	ref_file=$4

	files_1=$(find $source_dir/*_1.fastq.gz | head -${nr_samples})

        for file_1 in $files_1
        do
                base_file_1=$(basename "$file_1")
                # e.g. run=ERR525689

		#here we should check if the file already exists
		# "*.pair1.truncated.gz"

                run=${base_file_1%_1.fastq.gz}

		file_1_sample_dir=$sample_dir/${run}/${run}.pair1.truncated.gz

		if [ -f $file_1_sample_dir ]
		then
			echo $file_1_sample_dir "exists, no processing"
		else
			# echo $file_1_sample_dir "missing"
			run_sample $source_dir $sample_dir $run $ref_file
		fi
        done

	sample_stats $sample_dir
}

function run_all_samples_against_new_refs {

	sample_dir=$1
	nr_samples=$2
	ref_file=$3

	files_1=$(find $sample_dir/*/*pair1.truncated.gz | head -${nr_samples})

        for file_1 in $files_1
        do
                base_file_1=$(basename "$file_1")
                # e.g. run=ERR525689

                run=${base_file_1%.pair1.truncated.gz}
		#echo $run

		run_sample_against_new_refs $sample_dir $run $ref_file
	done

}

function sample_stats {
	paste	<(find $sample_dir/*/*idstats.txt |\
		while read LINE
		do
			base=$(basename "$LINE")
			sample=${base%.sorted.idstats.txt}
			echo $sample; \
		done) \
		<(cat $sample_dir/*/*idstats.txt | awk 'NR%2==1 {print $3}') \
		<(cat $sample_dir/*/*idstats.txt | awk 'NR%2==0 {print $4}') > $sample_dir/sample_stats.txt

	# cat $sample_dir/sample_stats.txt | sort -k2 -nr | awk '{print $0"\t"$2/$3}' > $sample_dir/sample_stats_sorted.txt

}

function run_sample {

	echo "start rum_sample"

        source_dir=$1
        sample_dir=$2
        run=$3
	ref_file=$4

	copy_and_run_trimming $source_dir $sample_dir $run

	refs=$(cat ${ref_file} | grep -v "#" | cut -f1 )

	for ref in $refs
	do
		echo "Starting to map sample against "$ref

		run_mapping $sample_dir $run $ref
		nr_mapped_reads=$(head -1 ${sample_dir}/${run}_${ref}/${run}.sorted.idstats.txt | cut -f3)

		if [ $nr_mapped_reads -gt 0 ]
		then
			echo "Nr of reads mapped: "$nr_mapped_reads
			run_diversiutils $sample_dir $run $ref

			run_calc_measures $sample_dir $run $ref
		else
			echo "No reads mapped to reference genome!"
		fi
	done
}


# this was my earlier sample to check for different results with untrimmed files: ERR1136746
# now I want to use ERR525804
function run_sample_against_new_refs {

	echo "start debug_sample"

        sample_dir=$1
        run=$2
	ref_file=$3

	# you do not have to copy and run trimming when running against new refs
	#copy_and_run_trimming $source_dir $sample_dir $run

	refs=$(cat ${ref_file} | grep -v "#" | cut -f1)

	for ref in $refs
	do
		echo "Starting to map sample against "$ref

		run_mapping $sample_dir $run $ref
		nr_mapped_reads=$(head -1 ${sample_dir}/${run}_${ref}/${run}.sorted.idstats.txt | cut -f3)

		if [ $nr_mapped_reads -gt 0 ]
		then
			echo "Nr of reads mapped: "$nr_mapped_reads
			run_diversiutils $sample_dir $run $ref

			run_calc_measures $sample_dir $run $ref
		else
			echo "No reads mapped to reference genome!"
		fi
	done
}

function copy_and_run_trimming {
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

	echo "start copying files from source dir"
	time /bin/cp -rfv $file_1 $sample_dir/${run}/
	time /bin/cp -rfv $file_2 $sample_dir/${run}/

	# redefine filenames, refer to sample_dir
        file_1=${sample_dir}/${run}/${run}_1.fastq.gz
        file_2=${sample_dir}/${run}/${run}_2.fastq.gz
	basename=${sample_dir}/${run}/${run}

	conda deactivate
	conda activate mgx

	echo "start adapter removal"
	time AdapterRemoval --file1 $file_1 --file2 $file_2 --basename $basename --trimns --trimqualities --collapse --threads 16 --minquality 25 --gzip

	rm -f $file_1
	rm -f $file_2

}

function run_mapping {

	#now make sure never to use the wrong version of bwa mem again
	conda deactivate
	conda activate mgx

	sample_dir=$1
	run=$2
	ref=$3

	#major fix: wrong second file!
	file_1=${sample_dir}/${run}/${run}.pair1.truncated.gz
	file_2=${sample_dir}/${run}/${run}.pair2.truncated.gz

	ll $file_1
	ll $file_2


	mkdir $sample_dir/${run}_${ref}
	#base for all other bam/same/etc files
	file=$sample_dir/${run}_${ref}/${run}

	ref_seq=scripts/mgx/ref_seqs/${ref}.fasta
	bwa index $ref_seq

        #TODO: check bwa mem options
	# bwa mem should work on .gz files
	echo "start bwa mem"

        time bwa mem -t 16 ${ref_seq} ${file_1} ${file_2} | samtools sort -@16 -o $file.sorted.bam -

	samtools view -b -F 4 $file.sorted.bam > $file.sorted.mapped.bam

	echo "start index sorted bam"
	#index bam
        time samtools index $file.sorted.bam
	time samtools index $file.sorted.mapped.bam

	echo "start idxstats"

        #create stats for mapped reads
        time samtools idxstats $file.sorted.bam > $file.sorted.idstats.txt
        time samtools idxstats $file.sorted.mapped.bam > $file.sorted.mapped.idstats.txt

	rm -f $file.sorted.bam
}

function run_diversiutils {

        sample_dir=$1
        run=$2
	ref=$3

	# to do: parametrize (and use naming convention to link both files)
	#ref=crassphage_refseq
	ref_seq=scripts/mgx/ref_seqs/${ref}.fasta
	coding_regions=scripts/mgx/ref_seqs/${ref}_codingregions.txt

	echo "start DiversiTools"

        time tools/DiversiTools/bin/diversiutils_linux -bam $sample_dir/${run}_${ref}/${run}.sorted.mapped.bam\
                -ref ${ref_seq}\
                -orfs ${coding_regions}\
                -stub $sample_dir/${run}_${ref}/${run}

	echo "start cleaning up DiversiTools results"

        #remove the "<NA>" strings
        time sed -e 's/<NA>//g' $sample_dir/${run}_${ref}/${run}_AA.txt > $sample_dir/${run}_${ref}/${run}_AA_clean.txt

	rm -f $sample_dir/${run}_${ref}/${run}_AA.txt
}

function run_calc_measures {


	echo "start calc_measures"

        sample_dir=$1
        run=$2
	ref=$3

        conda deactivate
        conda activate python37

	/bin/cp -rfv source/phages/codon_syn_non_syn_probabilities.txt ${sample_dir}/

        #create measures (based on ${run}_AA_clean.txt)
        time python source/phages/CalcDiversiMeasures.py -d $sample_dir -s $run -r $ref
}
