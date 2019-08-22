#-----------------------
# view_stats
#
# processing the results from [virsorter -> prodigal]
#-----------------------

#this is messy stuff, we have to clean it up
#purpose: getting the genome_id, contig_id, unique phage_name back from the generated individual protein names
#this is so messy because it is the result of the [virsorter -> prodigal] pipeline part
function translate_table {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	#cut off the category and put ; in front of genome_id
	#cat $gene_dir/IP_translation.txt | sed -e 's/-cat_[0-9]*_[0-9]*/;\0/' | sed -e 's/[0-9]*_[0-9]*_;/;\0/'

	#this one is for the prophages (ontaining genes)
	#cat $gene_dir/IP_translation.txt | sed -e 's/-cat_[0-9]*_[0-9]*/;\0/' | sed -e 's/_gene_[0-9]*_gene_[0-9]*-[0-9]*-[0-9]*;/;\0/g'

	#now do everything at the same time
	#we put in some ";" delimiters in intermediary steps to help selecting parts later on
	cat $gene_dir/IP_translation.txt | sed -e 's/-cat_[0-9]*_[0-9]*/;\0/' | sed -e 's/[0-9]*_[0-9]*_;/;\0/' |\
		sed -e 's/_gene_[0-9]*_gene_[0-9]*-[0-9]*-[0-9]*;/;\0/g' |\
		sed -e 's/___\([0-9]*_[0-9]*_\);/;\1/' |\
		sed -e 's/___\([0-9]*_[0-9]*_-circular\);/;\1/' |\
		sed -s 's/circular-cat/circular;-cat/g'

}

function genomes_from_proteins {
	translate_table | cut -d ';' -f2 | sed -e 's/_-circular$//g' | sed -e 's/__/_;_/g' | cut -d ';' -f1 | sed -e 's/_$//g'
}

function gene_part_from_proteins {
	translate_table | cut -d ';' -f2 | sed -e 's/\(gene_[0-9]*_gene_[0-9]*\)/;\1;/g' |\
		 sed -e 's/^[0-9]*_[0-9]*_//g' | sed -e 's/_;//' | cut -d ';' -f1
}

function contig_names_from_proteins {
	translate_table | cut -d ';' -f1 | cut -f2 | sed -e 's/VIRSorter_//g'
}

#we have the separate parts that make a phage unique and can optionally concatenate them again with _ (not strictly necessary)
function show_phage_names {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	#paste <(genomes_from_proteins) <(contig_names_from_proteins) <(gene_part_from_proteins) | awk '{print $1"_"$2"_"$3}'
	cat $gene_dir/IP_translation.txt | sed -e 's/_[0-9]*$/;\0/g' | cut -d ';' -f1 | cut -f2
}

function show_phages_with_most_proteins {
	show_phage_names | sort | uniq -c | sort -k1 -n -r
}

function show_genomes_with_most_phages {
	paste <(genomes_from_proteins) <(show_phage_names) | sort | uniq | cut -f1 | sort | uniq -c | sort -k1 -n -r
}

function show_phages_for_genome {
	genome=$1
	show_phage_names | grep $genome | sort | uniq
}


function stats_phages_1245 {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	echo " --- "
	echo "Nr of phages of categories 1,2,4,5 respectively:"
	show_phage_names | sort | uniq | grep -c cat_1
	show_phage_names | sort | uniq | grep -c cat_2
	show_phage_names | sort | uniq | grep -c cat_4
	show_phage_names | sort | uniq | grep -c cat_5
	echo " --- "
	echo "Nr of genes from phages of categories 1,2,4,5 respectively:"
	grep -c cat_1 $gene_dir/IP_translation.txt
	grep -c cat_2 $gene_dir/IP_translation.txt
	grep -c cat_4 $gene_dir/IP_translation.txt
	grep -c cat_5 $gene_dir/IP_translation.txt

}
