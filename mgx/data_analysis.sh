#ref_file=/hosts/linuxhome/mutant6/tmp/richard/ERP005989/ref_genome_ids.txt
#sample_dir=/hosts/linuxhome/mutant6/tmp/richard/ERP005989
#codon_table=source/phages/input_files/codon_syn_non_syn_probabilities.txt


function rerun_all_calc_measures {

	sample_dir=$1
	ref_file=$2

	refs=$(cat ${ref_file} | grep -v "#" | cut -f1)

        conda deactivate
        conda activate python37

	for ref in $refs
	do
		rerun_calc_measures_for_ref $sample_dir $ref
	done

	cat $sample_dir/*/*_sample_measures.txt | sort -r | uniq > $sample_dir/sample_measures.txt

}


function rerun_calc_measures_for_ref {

	sample_dir=$1
	ref=$2

        /bin/cp -rfv source/phages/codon_syn_non_syn_probabilities.txt ${sample_dir}/

        #create measures (based on ${run}_AA_clean.txt)
        time python source/phages/CalcDiversiMeasures.py -d $sample_dir -a -r $ref
}

function agg_sample_measures {

	sample_dir=$1
	# quick and dirty, get one header by reversing the order and uniq-ing the result
	cat $sample_dir/*/*_sample_measures.txt | sort -r | uniq > $sample_dir/sample_measures.txt
}


function run_all_gene_plots {

	sample_dir=$1
	ref_file=$2

	refs=$(cat ${ref_file} | grep -v "#" | cut -f1)

        conda deactivate
        conda activate python37

	for ref in $refs
	do
		echo $ref
		run_gene_plots_for_ref $sample_dir $ref 10 0.50
		run_gene_plots_for_ref $sample_dir $ref 10 0.95
	done
}


function run_gene_plots_for_ref {

	sample_dir=$1
	ref=$2
	depth=$3
	breadth=$4

	ref_seqs=scripts/mgx/ref_seqs

	conda deactivate
	conda activate python37

	# depth=10
	# breadth=0.20
	time python source/phages/MakeGenePlots.py -d $sample_dir -rd $ref_seqs -r $ref -td $depth -tb $breadth -ns 2
}

function run_family_gene_plots {

	sample_dir=$1
	breadth=$2

	ref_seqs=scripts/mgx/ref_seqs

	conda deactivate
	conda activate python37

	# depth=10
	# breadth=0.20
	# time python source/phages/MakeGenePlots.py -d $sample_dir -f -rd $ref_seqs -td 1 -tb $breadth
	time python source/phages/MakeGenePlots.py -d $sample_dir -f -rd $ref_seqs -td 10 -tb $breadth
}


function run_all_codon_measures {

	sample_dir=$1
	ref_file=$2
	codon_table=$3

	refs=$(cat ${ref_file} | grep -v "#" | cut -f1)

        conda deactivate
        conda activate python37

	for ref in $refs
	do
		python source/phages/CalcCodonMeasures.py -d $sample_dir -r $ref -c $codon_table
	done

}

function sample_stats {

	sample_dir=$1

	if [ -z "$sample_dir" ]
	then
		echo "please supply sample_dir as argument"
	else
		echo "writing sample stats to " ${sample_dir}/sample_stats.txt
		paste <(find $sample_dir/*/*sorted.idstats.txt |\
	                while read LINE
	                do
	                        base=$(basename "$LINE")
	                        sample=${base%.sorted.idstats.txt}
	                       echo $sample; \
        		done) \
		<(cat $sample_dir/*/*sorted.idstats.txt | awk 'NR%2==1 {print $1}') \
		<(cat $sample_dir/*/*sorted.idstats.txt | awk 'NR%2==1 {print $3}') \
		<(cat $sample_dir/*/*sorted.idstats.txt | awk 'NR%2==0 {print $4}') > $sample_dir/sample_stats.txt
	fi
}


function show_top_10 {
	sample_dir=$1
	echo "top 10 abundances of ref genomes in samples"
	echo "--"
	cat $sample_dir/sample_stats.txt | sort -k3 -nr | head -10
}

function show_all {
	sample_dir=$1
	cat $sample_dir/sample_stats.txt | sort -k3 -nr
}


function total_processed {

	sample_dir=/hosts/linuxhome/mutant6/tmp/richard/ERP005989

	ll ${sample_dir}/*/*.sorted.idstats.txt | wc -l
	sample_stats $sample_dir
	show_top_10 $sample_dir
}
