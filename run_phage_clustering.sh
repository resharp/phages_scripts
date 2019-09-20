#example: run_phage_clustering_workflow $gene_dir "0.5"
function run_phage_clustering_workflow {

	gene_dir=$1

	cut_off=$2

	#this one needs the output of PhageSharedContent.calc_shared_pc_measures()
	# e.g. shared_pc_measures.mcl_75.I25.txt
	# phages need at least 5 IPs
	filter_measures_on_nr_of_ips $gene_dir

	#abc is the simple label format for mcl (see https://micans.org/mcl/man/mcl.html#examples)
	#this also uses a cut-off of a minimum of 0.1 Jaccard distance
	make_abc_from_filtered_phage_measures $gene_dir $cut_off

	#run mcl with the Jaccard distance between 0.5 and 1 as weight
	#and inflation 2.5
	run_mcl_phages $gene_dir "2.5" $cut_off
	#run_mcl_phages $gene_dir "1.5 2.0 2.5 3.0"

	#make a table of { PHC_Id, PH_Id } coupling between phage clusters and phages
	split_mcl_phages $gene_dir "2.5" $cut_off

}

function show_phage_measure_stats {

	#gene_dir="../../_tools/mcl/1000"
	gene_dir=$1


	#TODO: remove grep -v of first line
	echo 'number of edges with jaccard index >= 0.5'

	nr_edges=$(cat ${gene_dir}/shared_pc_measures.mcl_75.I25.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.5) print $0}' | grep -v phage_1 | wc -l)
	nr_edges=$((nr_edges/2))	#prevent double counts
	echo $nr_edges

	echo 'number of nodes'
	cat ${gene_dir}/shared_pc_measures.mcl_75.I25.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.5) print $0}' | grep -v phage_1 | cut -d "," -f1 | sort | uniq | wc -l
}

function make_sif_from_phage_measures {

	#gene_dir="../../_tools/mcl/1000"
	gene_dir=$1

	out_name=shared_pc_measures.mcl_75.I25.sharing_0_1.sif

	cat ${gene_dir}/shared_pc_measures.mcl_75.I25.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.1) print $1,"sharing_0_1",$2}' | grep -v phage_1 > $gene_dir/$out_name

}

function filter_measures_on_nr_of_ips {
	gene_dir=$1

	out_name=shared_pc_measures.mcl_75.I25.nr_pc_min_5.txt
	cat ${gene_dir}/shared_pc_measures.mcl_75.I25.txt | awk ' BEGIN { FS=","} { if ($3 > 4 && $4 > 4) print $0 }' > $gene_dir/$out_name
}

function show_filtered_phage_measure_stats {

	gene_dir=$1


	echo 'number of edges with jaccard index >= 0.5 of phages containing at least 5 pcs'

	nr_edges=$(cat ${gene_dir}/shared_pc_measures.mcl_75.I25.nr_pc_min_5.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.5) print $0}' | wc -l)
	nr_edges=$((nr_edges/2))	#prevent double counts
	echo $nr_edges

	echo 'number of nodes (phages) containing at least 5 pcs with any jaccard index >= 0.5'
	cat ${gene_dir}/shared_pc_measures.mcl_75.I25.nr_pc_min_5.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.5) print $0}' | cut -d "," -f1 | sort | uniq | wc -l
}


#the sif can be used for visual inspection in Cytoscape
#this one uses a Jaccard index cut-off of 0.5
function make_sif_from_filtered_phage_measures {

	gene_dir=$1

	out_name=shared_pc_measures.mcl_75.I25.nr_pc_min_5.sharing_0_5.sif

	cat ${gene_dir}/shared_pc_measures.mcl_75.I25.nr_pc_min_5.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.5) print $1,"sharing_0_5",$2}' > $gene_dir/$out_name
}

#example: make_abc_from_filtered_phage_measures $gene_dir "0.5"
function make_abc_from_filtered_phage_measures {

	gene_dir=$1
	cut_off=$2

	cut_off_str=$(echo $cut_off | sed -e "s/\.//g")
	out_name=shared_pc_measures.mcl_75.I25.nr_pc_min_5.sharing_${cut_off_str}.abc

	#we filter on Jaccard index >= cut_off
	cat ${gene_dir}/shared_pc_measures.mcl_75.I25.nr_pc_min_5.txt | awk -v cut_off=$cut_off ' BEGIN { FS=","} { if ($6 >= cut_off) print $1,$2,$6}' > $gene_dir/$out_name
}

#example: view_abc_from_filtered_phage_measures $gene_dir "0.5"
function view_abc_from_filtered_phage_measures {
	gene_dir=$1
	cut_off=$2

	#we filter on Jaccard index >= cut_off
	cat ${gene_dir}/shared_pc_measures.mcl_75.I25.nr_pc_min_5.txt | awk -v cut_off=$cut_off ' BEGIN { FS=","} { if ($6 >= cut_off) print $1,$2,$6}'
}

#example:
# run_mcl_phages $gene_dir "1.5 2.0 2.5 3.0" "0.5"
function run_mcl_phages {

	gene_dir=$1

        #samples="1.5 2.0 2.5 3.0 5.0 8.0"
	samples=$2

	cut_off=$3

	cut_off_str=$(echo $cut_off | sed -e "s/\.//g")
	file_in_mcl=shared_pc_measures.mcl_75.I25.nr_pc_min_5.sharing_${cut_off_str}.abc

	#out_name=$gene_dir/out.shared_pc_measures.mcl_75.I25.nr_pc_min_5.sharing_0_1.abc

        for sample in $samples:
        do
                #added --expect-values for interpreting weight. Does it help?
                mcl $gene_dir/$file_in_mcl --abc -I ${sample} -te 4 --d --expect-values
        done

        for sample in $samples:
	do
		sample_str=$(echo "I"$sample | sed -e "s/\.//g")
		echo "largest clusters for inflation: "$sample_str
                cat $gene_dir/out.$file_in_mcl.$sample_str | awk '{print NF" phages"}' | head -3
        done
}

#example: split_mcl_phages $gene_dir "2.5" "0.5"
function split_mcl_phages {

        gene_dir=$1

        #samples denote runs for different mcl inflation factors
        sample=$2
	sample_str=$(echo "I"$sample | sed -e "s/\.//g")


	cut_off=$3
	cut_off_str=$(echo $cut_off | sed -e "s/\.//g")

        file_out_mcl=out.shared_pc_measures.mcl_75.I25.nr_pc_min_5.sharing_${cut_off_str}.abc.$sample_str
        file_out_split=split.shared_pc_measures.mcl_75.I25.nr_pc_min_5.sharing_${cut_off_str}.abc.$sample_str

        n=1
        cat $gene_dir/$file_out_mcl |\
                while read line
                do
                        for word in $line
                        do
                                echo PHC_$n $word >> $gene_dir/$file_out_split
                        done
                        n=$((n+1))
                done
}
