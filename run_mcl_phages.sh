#mcl is installed in conda env virsorter
function run_mcl_phages {

	#default inflation is 2.0. we run with 2.5 for comparing with cytoscape plug-in
	#keep it simple run with all default options

	#TODO: change location, this one is only running on 1000 genomes
	#---- * *  *     *               *
	gene_dir=/hosts/linuxhome/mgx/DB/PATRIC/patric/phage_genes_1000

	file_in_mcl=shared_gene_content.mcl_75.I20.txt

	#run with different inflations

	samples="1.5 2.0 2.5 3.0"
	#samples="2.0"
	for sample in $samples:
	do
		#added --expect-values for interpreting weight. Does it help?
		mcl $gene_dir/$file_in_mcl --abc -I ${sample} -te 4 --d --expect-values
	done

	samples="I15 I20 I25 I30"
	#samples="I20"

	for sample in $samples
		do echo "largest clusters for inflation: "$sample
		cat $gene_dir/out.$file_in_mcl.$sample | awk '{print NF" phages"}' | head -3
	done

}
