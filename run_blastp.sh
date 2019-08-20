
#make blast database
function make_blastp_db {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes

	#sample of some protein files
	#TODO: change to all proteins
	#----
	p_files=$(find $gene_dir/*/*.faa | head -1000)

	cat $p_files > $gene_dir/gene_samples.fasta

	#make translation file for translation of simple headers
	cat $gene_dir/gene_samples.fasta | sed -e 's/\s.*$//' |\
		awk 'BEGIN { num=1} { if($0~">") {print "IP_"num"\t"$0;num=num+1 }}' |\
		sed -e 's/>//g' > $gene_dir/IP_translation.txt

	#translate headers into simple >IP_[number] format
	cat $gene_dir/gene_samples.fasta | sed -e 's/\s.*$//' |\
		awk 'BEGIN { num=1} { if($0~">") {print ">IP_"num;num=num+1} else print $0  }' > $gene_dir/gene_samples_simple.fasta

	#makeblastdb -in gene_samples.fasta -parse_seqids -blastdb_version 5 -title "Phage gene samples" -dbtype prot out -out my.database_name
	makeblastdb -in $gene_dir/gene_samples_simple.fasta -parse_seqids -title "Phage gene samples" -dbtype prot -out $gene_dir/gene_samples_simple.faa

}

function cleanup_blastp_db {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes

	rm -f $gene_dir/*.fasta

	#remove fasta database
	rm -f $gene_dir/*.faa.*

	#remove pairwise blastp results
	rm -f $gene_dir/*.txt
}



#run blast against a fasta file
function run_blastp {

	evalue=$1

	#TODO: check if evalue is filled!

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes

	#blastp -db  -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

	echo 'start running blastp'
	blastp -query $gene_dir/gene_samples_simple.fasta -db $gene_dir/gene_samples_simple.faa\
		-evalue $evalue\
		-out $gene_dir/pairwise_blastout.txt\
		-num_threads 2\
		 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
	#-outfmt '7 qseqid sseqid length qlen qstart qend sstart send evalue'

	echo 'filtering out protein comparisons against themselves'
	cat $gene_dir/pairwise_blastout.txt | awk '{if ($1!=$2) print $0}' > $gene_dir/pairwise_blastout_filtered.txt

	echo 'number of hits:'
	wc -l $gene_dir/pairwise_blastout_filtered.txt

}


#view top 10 proteins within samples
function view_top_10_proteins {
	cat annotations/pairwise_blastout_filtered.txt | cut -f1 | sort | uniq -c | sort -k1 -n -r | head -10
}


function prepare_sif_for_cytoscape {

	#http://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html
	#sif format example:
	#node1 typeA node2
	#node2 typeB node3 node4 node5
	cat annotations/pairwise_blastout_filtered.txt | awk '{print $1 " homolog " $2}' > annotations/pw_blastout_cytoscape.sif


	#prepare format for mcl (without weight)
	cat annotations/pairwise_blastout_filtered.txt | awk '{print $1 "\t" $2}' > annotations/pw_blastout_cytoscape.abc

}

function run_mcl {

	#default inflation is 2.0. we run with 2.5 for comparing with cytoscape plug-in
	#keep it simple run with all default options
	#how 

	#run with different inflations

	samples="1.5 2.0 2.5 3.0"
	for sample in $samples:
	do
		mcl annotations/pw_blastout_cytoscape.abc --abc -I ${sample} -te 4 --d
	done

	samples="I15 I20 I25 I30"

	for sample in $samples
		do echo "largest clusters for inflation: "$sample
		cat annotations/out.pw_blastout_cytoscape.abc.$sample | awk '{print NF" proteins"}' | head -3
	done

}


#-----------------------
# TODO: move the following stuff, it has nothing to do with running blast, but processing the results from [virsorter -> prodigal]
#-----------------------

#this is messy stuff, we have to clean it up
#purpose: getting the genome_id, contig_id, unique phage_name back from the generated individual protein names
#this is so messy because it is the result of the [virsorter -> prodigal] pipeline part
function translate_table {

	#cut off the category and put ; in front of genome_id
	#cat annotations/IP_translation.txt | sed -e 's/-cat_[0-9]*_[0-9]*/;\0/' | sed -e 's/[0-9]*_[0-9]*_;/;\0/'

	#this one is for the prophages (ontaining genes)
	#cat annotations/IP_translation.txt | sed -e 's/-cat_[0-9]*_[0-9]*/;\0/' | sed -e 's/_gene_[0-9]*_gene_[0-9]*-[0-9]*-[0-9]*;/;\0/g'

	#now do everything at the same time
	#we put in some ";" delimiters in intermediary steps to help selecting parts later on
	cat annotations/IP_translation.txt | sed -e 's/-cat_[0-9]*_[0-9]*/;\0/' | sed -e 's/[0-9]*_[0-9]*_;/;\0/' |\
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
	#paste <(genomes_from_proteins) <(contig_names_from_proteins) <(gene_part_from_proteins) | awk '{print $1"_"$2"_"$3}'
	cat annotations/IP_translation.txt | sed -e 's/_[0-9]*$/;\0/g' | cut -d ';' -f1 | cut -f2
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
