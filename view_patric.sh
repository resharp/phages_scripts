#Some code to view and manipulate the genomes

function view_genome {

	genome=$1
	genome=$(echo $genome | tr '_' '.')

	less /hosts/linuxhome/mgx/DB/PATRIC/patric/ftp.patricbrc.org/genomes/$genome/$genome.fna

	#patric gene predictions are here:
	# /PATRIC.spgene.tab

}


#we will place all the phage predictions

# /hosts/linuxhome/mgx/DB/PATRIC/patric/phages

# and use the same structure
# /hosts/linuxhome/mgx/DB/PATRIC/patric/phages/$genome


#the files are now in

#mutant=mutant14
#/hosts/linuxhome/${mutant}/tmp/richard/virsorter_output/*/Predicted_viral_sequences/VIRSorter_cat-${category}.fasta
#/hosts/linuxhome/${mutant}/tmp/richard/virsorter_output/



function total_processed_genomes {

	total1=$(ll /hosts/linuxhome/mutant14/tmp/richard/virsorter_output/ | awk '{print $9}' | grep -E '[0-9]+\.[0-9]+' | wc -l)
	total2=$(ll /hosts/linuxhome/mutant31/tmp/richard/virsorter_output/ | awk '{print $9}' | grep -E '[0-9]+\.[0-9]+' | wc -l)
	total=$((total1+total2))
	echo $total
}


function copy_phage {

	mutant=$1
	genome=$2

	genome_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phages/$genome

	non_empty_files=$(ll /hosts/linuxhome/$mutant/tmp/richard/virsorter_output/$genome/Predicted_viral_sequences | grep cat | awk '{if($5!="0")print $9}')

	if [ ! -z "$non_empty_files" ]
	then
		mkdir $genome_dir
		for file in $non_empty_files
			#do echo /hosts/linuxhome/$mutant/tmp/richard/virsorter_output/$genome/Predicted_viral_sequences/$file
			do cp /hosts/linuxhome/$mutant/tmp/richard/virsorter_output/$genome/Predicted_viral_sequences/$file $genome_dir
		done
	fi
}


function copy_phages_from_mutant {

	mutant=$1
	genomes=$(ll /hosts/linuxhome/${mutant}/tmp/richard/virsorter_output/ | awk '{print $9}' | grep -E '[0-9]+\.[0-9]+')
	for genome in $genomes
		do copy_phage $mutant $genome
	done
}

