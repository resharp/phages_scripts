source scripts/job_stats.sh


function get_phage_filenames_for_cat {

	mutant=$1
	category=$2

        if [[ "$category" =~ ^("1"|"2"|"3")$ ]]
        then
		ll /hosts/linuxhome/${mutant}/tmp/richard/virsorter_output/*/Predicted_viral_sequences/VIRSorter_cat-${category}.fasta |\
			awk '{if($5!="0")print $9}'
        else
                ll /hosts/linuxhome/${mutant}/tmp/richard/virsorter_output/*/Predicted_viral_sequences/VIRSorter_prophages_cat-${category}.fasta |\
                        awk '{if($5!="0")print $9}'
        fi

	#TODO: filter out non-zero files

	#here are all samples (not just the ones containing phages)
	# ll /hosts/linuxhome/${mutant}/tmp/richard/virsorter_output/*/Predicted_viral_sequences/VIRSorter_cat-${category}.fasta | awk '{if($5!="0")print $9}'

	#samples=$(ll /hosts/linuxhome/mutant31/tmp/richard/virsorter_output/ | tr ' ' '|' |  grep -v "^total" | grep 4096 | grep -v "\.\." | cut -d "|" -f15-)

	#now we probably want to join with the non-zero files

	#echo $samples
}

function get_phage_samples_for_cat {

	mutant=$1
	category=$2

        if [[ "$category" =~ ^("1"|"2"|"3")$ ]]
        then
		ll /hosts/linuxhome/${mutant}/tmp/richard/virsorter_output/*/Predicted_viral_sequences/VIRSorter_cat-${category}.fasta |\
			awk '{if($5!="0")print $9}' |\
        		cut -d "/" -f8
	else
                ll /hosts/linuxhome/${mutant}/tmp/richard/virsorter_output/*/Predicted_viral_sequences/VIRSorter_prophages_cat-${category}.fasta |\
                        awk '{if($5!="0")print $9}' |\
			cut -d "/" -f8
        fi
}

#temporary function
function get_fasta_file {
	cat_fasta mutant14 1 | head -1 | tr ' ' '|' | cut -d '|' -f11
}

function run_prodigal {

	fasta_file=$(get_fasta_file)

	echo "Run prodigal on" $fasta_file

	prodigal -i $fasta_file -o annotations/my.genes -a annotations/my.proteins.faa -p meta
}

function run_prodigal_on_some_samples {

	mutant=mutant14;samples=$(get_phage_samples_for_cat $mutant 1 | head -10)

	#only category 1 and 2 (and 4 and 5)
	for sample in $samples
		do prodigal -i <(cat /hosts/linuxhome/$mutant/tmp/richard/virsorter_output/1797198.3/Predicted_viral_sequences/VIRSorter_*-[1-2\|4-5].fasta)\
			-o annotations/$sample.genes -a annotations/$sample.proteins.faa -p meta
	done
}

