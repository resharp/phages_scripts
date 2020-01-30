 # run complete pipeline on one SRA project

# e.g. run_all_samples $source_dir $sample_dir 10
source_dir=/hosts/linuxhome/mutant14/tmp/richard/sra/ERP005989
sample_dir=/hosts/linuxhome/mutant23/tmp/richard/ERP005989

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

	# cp $file_1 $sample_dir/${basename}_1.fastq.gz
	/bin/cp -rfv $file_1 $sample_dir/
	/bin/cp -rfv $file_2 $sample_dir/

	# redefine filenames, refer to sample_dir
        file_1=${sample_dir}/${run}_1.fastq.gz
        file_2=${sample_dir}/${run}_2.fastq.gz
	basename=${sample_dir}/${run}

	conda deactivate
	conda activate mgx

	AdapterRemoval --file1 $file_1 --file2 $file_2 --basename $basename --trimns --trimqualities --collapse --threads 16 --minquality 25 --gzip

}

function run_mapping {

	sample_dir=$1
	run=$2
	file_1=${sample_dir}/${run}.pair1.truncated.gz
	file_2=${sample_dir}/${run}.pair1.truncated.gz

	ll $file_1
	ll $file_2

	#base for all other bam/same/etc files
	mkdir $sample_dir/${run}
	file=$sample_dir/${run}/${run}

	ref_seq=scripts/mgx/ref_seqs/crassphage_refseq.fasta

        #TODO: check bwa mem options
	# bwa mem should work on .gz files
        bwa mem -t 16 ${ref_seq} ${file_1} ${file_2} > ${file}.sam

	#bwa mem -t 24 /linuxhome/tmp/stijn/SILVA_HLA/SILVA_HLA.fasta /linuxhome/tmp/stijn/SRA_trimmed/PE/${base}.pair1.truncated.gz \
	#	/linuxhome/tmp/stijn/SRA_trimmed/PE/${base}.pair2.truncated.gz

	#TODO: check if you can pipe output directly to prevent space usage
        #convert sam to bam
        samtools view -@ 16 -S -b $file.sam > $file.bam

        #sort bam
        samtools sort -@ 16 -o $file.sorted.bam $file.bam

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
