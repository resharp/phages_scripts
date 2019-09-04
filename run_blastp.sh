#example workflow
#depends on prodigal results (see run_prodigal.sh)
function run_blastp_workflow {
	
	#TODO: split two kind of working directories
	# one where the *.faa files are (original gene predictions from prodigal)
	# and one to further process results
	
	#also builds the IP_translation.txt table that is heavily used in extracting relations in view_stats.sh
	make_blastp_db
	
	run_blastp 1e-3

	filter_pairwise_blastp_hits_on_coverage

}


#make blast database
function make_blastp_db {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes

	#concat all files!
	p_files=$(find $gene_dir/*/*.faa)

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

#this function moves the results of a subset of 1000 genomes to a new directory for analysis
function move_blastp_results {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes
	to_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	mv -f $gene_dir/*.fasta $to_dir

	#remove fasta database
	mv -f $gene_dir/*.faa.* $to_dir

	#remove pairwise blastp results
	mv -f $gene_dir/*.txt $to_dir
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

function filter_pairwise_blastp_hits_on_coverage {

	#TODO: change to real gene_dir, not the 1000-sample
	# * *  *   *    *
	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_all

	#now we want to filter on 75% pairwise query and 75% target coverage
	#calculate (to-from+1)/(length) both for query and for target
	cat $gene_dir/pairwise_blastout_filtered.txt |\
		 awk '{ print $0"\t"($8-$7+1)/$13"\t"($10-$9+1)/$14}' |\
		 awk '{if($15 > 0.75 && $16 > 0.75) print $0}' > $gene_dir/pairwise_blastout_filtered_75coverage.txt
}


#view top 10 proteins within samples
function view_top_10_proteins {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes

	#cat $gene_dir/pairwise_blastout_filtered.txt | cut -f1 | sort | uniq -c | sort -k1 -n -r | head -10
	cat $gene_dir/pairwise_blastout.txt | cut -f1 | sort | uniq -c | sort -k1 -n -r | head -10
}

function view_average_nr_blastp_hits {

	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes

	cat $gene_dir/pairwise_blastout.txt | cut -f1 | sort | uniq -c | sed -e 's/ IP_[0-9]*$/;\0/g' |\
		 cut -d ';' -f1 | awk '{ total += $1; count++ } END { print total/count }'
}

function view_protein {
	protein=$1

	#grep -E "$protein\b" -A100 patric/patric/phage_genes/gene_samples_simple.fasta | less

	#seqtk is a module installed in the python37 conda env
	seqtk subseq $gene_dir/gene_samples_simple.fasta <(echo $protein) | less
}
