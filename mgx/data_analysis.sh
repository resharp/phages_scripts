#ref_file=/hosts/linuxhome/mutant6/tmp/richard/ERP005989/ref_genome_ids.txt
#sample_dir=/hosts/linuxhome/mutant6/tmp/richard/ERP005989

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
		run_gene_plots_for_ref $sample_dir $ref 10 0.20
		run_gene_plots_for_ref $sample_dir $ref 10 0.80
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
	time python source/phages/MakeGenePlots.py -d $sample_dir -f -rd $ref_seqs -td 1 -tb $breadth
	time python source/phages/MakeGenePlots.py -d $sample_dir -f -rd $ref_seqs -td 10 -tb $breadth
}
