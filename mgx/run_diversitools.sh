#DiversiTools
#http://josephhughes.github.io/DiversiTools/tutorial.html

sample=MGXDB008660
sample_dir=/hosts/linuxhome/mutant31/tmp/richard/crassphage_samples/${sample}

function run_diversiutils {

	sample=$1

	tools/DiversiTools/bin/diversiutils_linux -bam $sample_dir/${sample}_filtered.sorted.bam\
		-ref $sample_dir/crassphage_refseq.fasta\
		-orfs $sample_dir/codingregions_new.txt\
		-stub $sample_dir/${sample}
}

#cut -f3,5,11,12,13 $sample_dir/MGXDB008660_CDS_AA.txt | less

# $ perl bin/diversiutils.pl -bam path/to/input.bam -ref path/to/reference.fasta -orfs CodingRegions.txt -stub out

#Protein	Beg	End	Reference
#KP06_gp66	61171	62631	NC_024711.1
#KP06_gp74	83140	84612	NC_024711.1


#KP06_gp74
#CDS complement(83140..84612)
#/product="putative Major capsid protein"

#CDS complement(61171..62631)
#/product="putative Phage tail-collar fibre protein

function convert_aa_file {

	sample=$1

	#remove the "<NA>" strings
	sed -e 's/<NA>//g' $sample_dir/${sample}_AA.txt > $sample_dir/${sample}_AA_clean.txt

}


gene=KP06_gp66
#view_gene_metrics $gene
function view_gene_metrics {
	gene=$1
	# | less
}

