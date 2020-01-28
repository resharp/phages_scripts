
# AdapterRemoval on all paired files in a certain directory?
# the directory is

project_dir=/hosts/linuxhome/mutant14/tmp/richard/sra/ERP005989

# for each file in project_dir
function run_adapter_removal_on_project {
	project_dir=$1

	# ERR525689

	# to do: remove head -1 if you want to run all files
	files_1=$(find $project_dir/*_1.fastq.gz | head -1)

	for file_1 in $files_1
	do
		file_2=$(echo $file_1 | sed -e 's/_1/_2/g')
		run_adapter_removal $project_dir $file_1 $file_2
	done

}

function run_adapter_removal {

	project_dir=$1
	file_1=$2
	file_2=$3

	echo "Starting to run AdapterRemoval"
	echo $file_1
	echo $file_2
	base=ERR525689.output_paired

	# to do use --gzip option for smaller files
	# --gzip
	AdapterRemoval --file1 $file_1 --file2 $file_2 --basename $project_dir/$base --trimns --trimqualities --collapse --threads 2 --minquality 25 --gzip

}
