#-----------------------
# view_stats
#
# processing the results from [virsorter -> prodigal]
# but depends on IP_translation built in run_blastp.sh
#
# building tables for coupling of phage to IP (this could be moved to end of run_prodigal script)
#
# also translating gene content into protein clusters of a specific clustering
#	(this could be moved to end of run_mcl script)
#-----------------------

#this is messy stuff, we have to clean it up
#purpose: getting the genome_id, contig_id, unique phage_name back from the generated individual protein names
#this is so messy because it is the result of the [virsorter -> prodigal] pipeline part
function translate_table {

	gene_dir=$1

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
	gene_dir=$1
	translate_table $gene_dir | cut -d ';' -f2 | sed -e 's/_-circular$//g' | sed -e 's/__/_;_/g' | cut -d ';' -f1 | sed -e 's/_$//g'
}

function gene_part_from_proteins {
	gene_dir=$1
	translate_table $gene_dir | cut -d ';' -f2 | sed -e 's/\(gene_[0-9]*_gene_[0-9]*\)/;\1;/g' |\
		 sed -e 's/^[0-9]*_[0-9]*_//g' | sed -e 's/_;//' | cut -d ';' -f1
}

function contig_names_from_proteins {
	gene_dir=$1
	translate_table $gene_dir | cut -d ';' -f1 | cut -f2 | sed -e 's/VIRSorter_//g'
}

#example:
# show_phage_names /hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000
function show_phage_names_and_ips {

	gene_dir=$1

	#we have the separate parts that make a phage unique and can optionally concatenate them again with _ (not strictly necessary)
	#paste <(genomes_from_proteins) <(contig_names_from_proteins) <(gene_part_from_proteins) | awk '{print $1"_"$2"_"$3}'

	cat $gene_dir/IP_translation.txt | sed -e 's/_[0-9]*$/;\0/g' | cut -d ';' -f1 | cut -f1,2
}

function show_phage_names {

	gene_dir=$1

	#we have the separate parts that make a phage unique and can optionally concatenate them again with _ (not strictly necessary)
	#paste <(genomes_from_proteins) <(contig_names_from_proteins) <(gene_part_from_proteins) | awk '{print $1"_"$2"_"$3}'

	cat $gene_dir/IP_translation.txt | sed -e 's/_[0-9]*$/;\0/g' | cut -d ';' -f1 | cut -f2
}


# format: PH_Id, PH_Name
# from prodigal gene predictions (not specific to any clustering)
# depends on IP_translation.txt built in run_blastp.sh
function show_phage_table {

	gene_dir=$1

	#we have the separate parts that make a phage unique and can optionally concatenate them again with _ (not strictly necessary)
	#paste <(genomes_from_proteins) <(contig_names_from_proteins) <(gene_part_from_proteins) | awk '{print $1"_"$2"_"$3}'

	#sort on phage name and add a unique PH_Id number
	cat $gene_dir/IP_translation.txt | sed -e 's/_[0-9]*$/;\0/g' | cut -d ';' -f1 |\
		cut -f2 | sort | uniq |\
	 	awk 'BEGIN { num=1} { print "PH_"num"\t"$0;num=num+1 }'
}


# format: PH_Id, IP_Id, PH_Name
function make_phage_ip_table {

	gene_dir=$1

	#join on phage_name (show_phage_table already sorted on phage name)
	join -1 2 -2 2 -o2.1,1.1,0 <(cat $gene_dir/IP_translation.txt | sed -e 's/_[0-9]*$/;\0/g' | cut -d ';' -f1 | sort -k2) <(show_phage_table $gene_dir)\
		 > $gene_dir/phage_ip_table.txt
}

function make_phage_ip_table_short {

	gene_dir=$1

	cat $gene_dir/phage_ip_table.txt | cut -d " " -f1,2 > $gene_dir/phage_ip_table_short.txt
}


#format: PH_Id, PH_Name
#depends on phage_ip_table
function make_phage_table {

	gene_dir=$1

	cat $gene_dir/phage_ip_table.txt | awk '{ print $1"\t"$3}' | uniq > $gene_dir/phage_table.txt
}


# This functions uses a specific clustering result
# e.g. mcl_75.I20 (this means 75% query and target coverage threshold and inflation factor 2.0)
#format: IP_Id, PC_Id
#example: make_ip_pc_table $gene_dir I20
function make_ip_pc_table {

	gene_dir=$1
	
	sample=$2

	pc_file=out.pw_blastout_mcl_75.abc.$sample

	ip_pc_table=$gene_dir/ip_pc_table.mcl_75.$sample
	pc_table=$gene_dir/pc_table.mcl_75.$sample

	#take for example
	# out.pw_blastout_mcl_75.abc.I20

	#then put PC_[nr] in front of each line

	#then loop through lines
	#loop through words
	#first word is PC_Id
	#for every other word
	#write IP_Id, PC_Id

	>$pc_table
	cat $gene_dir/$pc_file | awk 'BEGIN { num=1} { print "PC_"num;num=num+1 }' > $pc_table

	> $ip_pc_table
	cat $gene_dir/$pc_file | awk 'BEGIN { num=1} { print "PC_"num"\t"$0;num=num+1 }' |\
 	        while read line
                do
			first=false
			for word in $line
			do
        	                if [ "$first" = "false" ]
				then
					pc=$word
					first=true
				else
					echo $word $pc >> $ip_pc_table
				fi
                	 done
		done

}

#now we make a smaller ip_pc table for a specific clustering for pc that contain at least 100 ips
#this is because we want to have a group of at least 100 phages, so e.g. 100 ips should be found at least for that pc
#ip=individual protein
#pc=protein cluster
#example
#make_filtered_ip_pc_table $gene_dir "mcl_75.I25" 100
function make_filtered_ip_pc_table {

	gene_dir=$1
	extension=$2
	min_nr_phages=$3

	out_file=$gene_dir/ip_pc_table.$extension.filter_$min_nr_phages

	#first determine the pcs that result from at least min_nr_phages ips
	nr_pcs=$(cat $gene_dir/ip_pc_table.$extension | cut -d " " -f2 | uniq -c |\
		 awk -v min_nr_phages=$min_nr_phages '{if ($1 > min_nr_phages-1 ) print $0}' | wc -l)

	echo "Nr of PCs :"$nr_pcs

	> $out_file
	cat $gene_dir/ip_pc_table.$extension |\
		awk -F_ '{print $3"\t"$0}'| awk -v nr_pcs=$nr_pcs '{ if ($1 < nr_pcs + 1 ) print $2"\t"$3}' > $out_file

}

function show_phages_with_most_proteins {
	gene_dir=$1

	cat $gene_dir/phage_ip_table.txt | tr " " "|" | cut -d '|' -f1 | uniq -c | sort -k1 -n -r

	#show_phage_names_and_ips $gene_dir | sort | uniq -c | sort -k1 -n -r
}

function show_genomes_with_most_phages {
	gene_dir=$1
	paste <(genomes_from_proteins $gene_dir ) <(show_phage_names $gene_dir) | sort | uniq | cut -f1 | sort | uniq -c | sort -k1 -n -r
}

function show_phages_for_genome {
	genome=$1
	gene_dir=$1
	show_phage_names $gene_dir | grep $genome | sort | uniq
}


#example
# stats_phages_1245 /hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000
function stats_phages_1245 {

	gene_dir=$1

	echo " --- "
	echo "Nr of phages of categories 1,2,4,5 respectively:"
	show_phage_names $gene_dir | sort | uniq | grep -c cat_1
	show_phage_names $gene_dir | sort | uniq | grep -c cat_2
	show_phage_names $gene_dir | sort | uniq | grep -c cat_4
	show_phage_names $gene_dir | sort | uniq | grep -c cat_5
	echo " --- "
	echo "Nr of genes from phages of categories 1,2,4,5 respectively:"
	grep -c cat_1 $gene_dir/IP_translation.txt
	grep -c cat_2 $gene_dir/IP_translation.txt
	grep -c cat_4 $gene_dir/IP_translation.txt
	grep -c cat_5 $gene_dir/IP_translation.txt

}
