#ref_file=/hosts/linuxhome/mutant6/tmp/richard/ERP005989/ref_genome_ids.txt
#sample_dir=/hosts/linuxhome/mutant6/tmp/richard/ERP005989

function rerun_all_calc_measures {

	sample_dir=$1
	ref_file=$2

	refs=$(cat ${ref_file} | grep -v "#")

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


