
#make blast database
function make_blastp_db {

	#sample of some protein files
	p_files=$(find patric/patric/phage_genes/*.faa | head -100)

	cat $p_files > annotations/gene_samples.fasta

	#make translation file for translation of simple headers
	cat annotations/gene_samples.fasta | sed -e 's/\s.*$//' |\
		awk 'BEGIN { num=1} { if($0~">") {print "IP_"num"\t"$0;num=num+1 }}' |\
		sed -e 's/>//g' > annotations/IP_translation.txt

	#translate headers into simple >IP_[number] format
	cat annotations/gene_samples.fasta | sed -e 's/\s.*$//' |\
		awk 'BEGIN { num=1} { if($0~">") {print ">IP_"num;num=num+1} else print $0  }' > annotations/gene_samples_simple.fasta

	#makeblastdb -in gene_samples.fasta -parse_seqids -blastdb_version 5 -title "Phage gene samples" -dbtype prot out -out my.database_name
	makeblastdb -in annotations/gene_samples_simple.fasta -parse_seqids -title "Phage gene samples" -dbtype prot -out annotations/gene_samples_simple.faa

}


#run blast against a fasta file
function run_blastp {

	evalue=$1

	#blastp -db  -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

	echo 'start running blastp'
	blastp -query annotations/gene_samples_simple.fasta -db annotations/gene_samples_simple.faa\
		-evalue $evalue\
		-out annotations/pairwise_blastout.txt\
		-num_threads 2\
		 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
	#-outfmt '7 qseqid sseqid length qlen qstart qend sstart send evalue'

	echo 'filtering out protein comparisons against themselves'
	cat annotations/pairwise_blastout.txt | awk '{if ($1!=$2) print $0}' > annotations/pairwise_blastout_filtered.txt

	echo 'number of hits:'
	wc -l annotations/pairwise_blastout_filtered.txt

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

