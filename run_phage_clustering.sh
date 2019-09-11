
function show_phage_measure_stats {
	gene_dir="../../_tools/mcl/1000"

	echo 'number of edges with jaccard index >= 0.5'
	
	nr_edges=$(cat ${gene_dir}/shared_pc_measures.mcl_75.I20.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.5) print $0}' | grep -v phage_1 | wc -l)
	nr_edges=$((nr_edges/2))	#prevent double counts
	echo $nr_edges
	
	echo 'number of nodes'
	cat ${gene_dir}/shared_pc_measures.mcl_75.I20.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.5) print $0}' | grep -v phage_1 | cut -d "," -f1 | sort | uniq | wc -l
}

function make_sif_from_phage_measures {

	gene_dir="../../_tools/mcl/1000"
	
	out_name=shared_pc_measures.mcl_75.I20.sharing_0_5.sif
	
	cat ${gene_dir}/shared_pc_measures.mcl_75.I20.txt | awk ' BEGIN { FS=","} { if ($6 >= 0.5) print $1,"sharing_0_5",$2}' | grep -v phage_1 > $gene_dir/$out_name
	
}
