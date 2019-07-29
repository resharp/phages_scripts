##search anything in all of the scripts in the scripts directory
function vs {

	target=$1
	lines_after=${2:-0}

	cat scripts/*.sh | grep $target -A ${lines_after} | less
}

function vh {

	target=$1
	lines_after=${2:-0}

	history | grep $target -A ${lines_after} | less
}


####################################################################
# functions for the single genome parallel runs
#
####################################################################

# show non empty signal files
# if all categories are 0 the size of the file is 1823. we don't want those
function non_empty_signal_files {
	mutant=$1

	ll /hosts/linuxhome/$mutant/tmp/richard/virsorter_output/[1-2]*/VIR*signal.csv | sort -k7 | awk '{if($5!="1823" && $5!="0")print $0}'

}

function cat_fasta {

	mutant=$1
	category=$2

	if [[ "$category" =~ ^("1"|"2"|"3")$ ]]
	then
		ll /hosts/linuxhome/$mutant/tmp/richard/virsorter_output/*/Predicted_viral_sequences/VIRSorter_cat-${category}.fasta |\
			awk '{if($5!="0")print $0}'
	else
		ll /hosts/linuxhome/$mutant/tmp/richard/virsorter_output/*/Predicted_viral_sequences/VIRSorter_prophages_cat-${category}.fasta |\
			awk '{if($5!="0")print $0}'
	fi
}

####################################################################
# functions for the runs on the concatenated files of 1000 genomes
#
####################################################################


#view ready
function vr {
	view_ready ~/virsorter_admin/virsorter_scheduler.txt
}

function view_size {
	# printf "My name is \"%s\".\nIt's a pleasure to meet "%s"\n" "John" "Mary"
	size=$(du -h $1 | cut -f1 | tail -1)

	printf %s" size of : "%s"\n" $size $1
}


#view the size of all output directories on the different mutants
function view_all_vs_out {

	rm -f view_all_vs_out.txt
	outputs=$(ll /hosts/linuxhome/mutant*/tmp/richard/virsorter_output/genome_ids_split_* | grep _split_ | sed -e 's/://g')
	for dir in $outputs; do view_size $dir >> view_all_vs_out.txt; done
	echo "------------"
	echo "total size: "
	cat view_all_vs_out.txt | tr ' ' '\t' | cut -f1 | sed -e 's/G//g' | awk 'BEGIN {} ; {sum+=$1} END {print sum}'

}

#view last ready of the 1000-batches
function vlr {
	cat ~/virsorter_admin/virsorter_scheduler.txt | grep READY | cut -f1,3 | awk '{ print "/hosts/linuxhome/"$2"/tmp/richard/virsorter_output/"$1 }'| \
	tr "\." "_" | cut -f2 | while read line; do ll "$line/VIRSorter_global-phage-signal.csv"; done | sort -k 7,8
}

#prepare statements for analyzing of the iterations
function prep_ana_iter {
	cat virsorter_admin/virsorter_scheduler.txt | grep READY | cut -f1,3 | \
	awk '{ print $1"\t"" ana_iterations /hosts/linuxhome/"$2"/tmp/richard/virsorter_output/"$1 " >> results/ana_iter_"$1}'| tr "\." "_" | cut -f2 | \
	while read line; do echo "$line"; done
}

function prep_show_vs {
	cat virsorter_admin/virsorter_scheduler.txt | grep READY | cut -f1,3 | \
	awk '{ print $1"\t"" show_vs /hosts/linuxhome/"$2"/tmp/richard/virsorter_output/"$1 " >> results/show_vs_"$1}'| tr "\." "_" | cut -f2 | \
	while read line; do echo "$line"; done
}

