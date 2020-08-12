# obsolete code, only used for pilot!
#
# DiversiTools
# http://josephhughes.github.io/DiversiTools/tutorial.html
#
# CalcDiversiMeasures.py
# https://github.com/resharp/phages
#
#

#sample=MGXDB008660

#TODO: change to other mutant
#------ ***
sample_dir=/hosts/linuxhome/mutant14/tmp/richard/crassphage_samples

master_dir=/hosts/linuxhome/chaperone/tmp/richard/crassphage_samples
admin_dir=/hosts/linuxhome/mgx/DB/MGXDB

function run_gene_samples {

	# run 20 gene samples
	#TODO: change to the number of samples you would like to use
	# samples=$(cat $admin_dir/taxon_1211417_counts_sorted.txt | awk '{if ($3>10000) print $1}' | head -20 | grep -v "sample")


	samples=$(grep NC_024711.1 $master_dir/MGXDB*/*.sorted.idstats.txt | sort -k 3 -n -r | awk '{if ($3 > 106000 && $3 < 170000) print $0}' | cut -f1 | cut -d "/" -f8)

	#samples=$(echo "MGXDB011779")

	for sample in $samples
	do
		echo $sample

		#copy bam file to $sample_dir

		mkdir $sample_dir/${sample}
		/bin/cp $master_dir/${sample}/${sample}_filtered.sorted.bam $sample_dir/${sample}/${sample}_filtered.sorted.bam
		/bin/cp $master_dir/${sample}/${sample}_filtered.sorted.bam.bai $sample_dir/${sample}/${sample}_filtered.sorted.bam.bai
		#ll $sample_dir/${sample}

		echo "start run diversiutils"
		run_diversiutils $sample_dir $sample
		#ll $sample_dir/${sample}

		run_calc_measures $sample_dir $sample

		#copy results back to right directory in master dir
		/bin/cp $sample_dir/${sample}/${sample}_*.txt $master_dir/${sample}


		#clean up (not all)
		rm -f $sample_dir/${sample}/${sample}_filtered.sorted.bam
		rm -f $sample_dir/${sample}/${sample}_filtered.sorted.bam.bai
		rm -f $sample_dir/${sample}/${sample}_AA.txt
		rm -f $sample_dir/${sample}/${sample}_entropy.txt
		rm -f $sample_dir/${sample}/${sample}_read.txt

	done
}


function run_diversiutils {

	sample_dir=$1
	sample=$2

	tools/DiversiTools/bin/diversiutils_linux -bam $sample_dir/${sample}/${sample}_filtered.sorted.bam\
		-ref $sample_dir/crassphage_refseq.fasta\
		-orfs scripts/mgx/crassphage_codingregions.txt\
		-stub $sample_dir/${sample}/${sample}

	#remove the "<NA>" strings
	sed -e 's/<NA>//g' $sample_dir/${sample}/${sample}_AA.txt > $sample_dir/${sample}/${sample}_AA_clean.txt
}

function run_calc_measures {

	sample_dir=$1
	sample=$2

	conda deactivate
	conda activate python37

	#create measures (based on ${sample}_AA_clean.txt)
	python source/phages/CalcDiversiMeasures.py -d $sample_dir -s $sample
}


function rerun_all_calc_measures {

	sample_dir = $1
	samples=$(ls $sample_dir | grep MGXDB)
	for sample in $samples
	do
		run_calc_measures $sample_dir $sample
	done

	python source/phages/MakeGenePlots.py -d $sample_dir
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
	sed -e 's/<NA>//g' $sample_dir/${sample}/${sample}_AA.txt > $sample_dir/${sample}/${sample}_AA_clean.txt

}

