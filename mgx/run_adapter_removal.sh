
# AdapterRemoval on all paired files in a certain directory
# the directory is

project_dir=/hosts/linuxhome/mutant14/tmp/richard/sra/ERP005989

# for each file in project_dir
# Adapter removal and quality filtering of the paired-end fastq files was done with AdapterRemoval v2
function run_adapter_removal_on_project {
	project_dir=$1

	# ERR525689

	# to do: remove head -1 if you want to run all files
	files_1=$(find $project_dir/*_1.fastq.gz | head -4)

	for file_1 in $files_1
	do
		base_file_1=$(basename "$file_1")
		# e.g. run=ERR525689

		run=${base_file_1%_1.fastq.gz}

		run_adapter_removal $project_dir $run
	done

}

function run_adapter_removal {

	project_dir=$1
	run=$2

	file_1=${project_dir}/${run}_1.fastq.gz
	file_2=${project_dir}/${run}_2.fastq.gz
	basename=${project_dir}/${run}

	echo "Starting to run AdapterRemoval for "$run" in "$project_dir
	#echo $file_1
	#echo $file_2
	#echo $basename
	# e.g. basename=${project_dir}/ERR525689

	AdapterRemoval --file1 $file_1 --file2 $file_2 --basename $basename --trimns --trimqualities --collapse --threads 16 --minquality 25 --gzip

}
