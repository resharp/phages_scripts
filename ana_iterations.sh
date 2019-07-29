function analyze_one_iteration {

	out_dir=$1
	iteration=$2
	nr_of_genes_to_be_shown=$3

	echo "Start analyzing $iteration ..."
	echo "------------------------------"

	##The result of the first iteration is in r_1
	##We want to determine the pangenome. Is it large enough for clustering?

	##determine sequences of category 1 or 2 (583)
	## Let us just store the sequence names in /dev/shm/phages_1_or_2.txt
	rm -f /dev/shm/phages_1_or_2.txt
	cat $out_dir/$iteration/phage_signal_backup.csv | awk '{ if(( $6 == 1) || ($6 == 2)) print $1}' >> /dev/shm/phages_1_or_2.txt
	## less /dev/shm/phages_1_or_2.txt

	##Build all combinations of phages and gene predictions in the order of phages and then genes
	rm -f /dev/shm/all_phages_plus_phage_genes.txt
	cat $out_dir/$iteration/affi_backup.csv | grep -v '>' | tr '|' '\t' | awk '{print $1,"\t",$6}' | awk '{if($2!="-" && $2!="")print $0}' | sed -e 's/\-gene_[0-9]*//g' | sort -k 1,1 -k2,2 | uniq >> /dev/shm/all_phages_plus_phage_genes.txt
	## less /dev/shm/all_phages_plus_phage_genes.txt

	##now intersect with the cat1,2 phages (this step is unnecessary for cat 1 and 2, because cat 1 and 2 always contain gene predictions)
	rm -f /dev/shm/phages_that_contain_phage_genes.txt
	join -1 1 -2 1 <(cat /dev/shm/phages_1_or_2.txt| sort)  <(cat /dev/shm/all_phages_plus_phage_genes.txt | cut -f1 | sort | uniq) \
		>> /dev/shm/phages_that_contain_phage_genes.txt
	## less /dev/shm/phages_that_contain_phage_genes.txt
	echo "Nr of phages of category 1 or 2:"
	wc -l /dev/shm/phages_that_contain_phage_genes.txt

	##now we only want the phages that have at least one gene in the shared_phage_genes
	rm -f /dev/shm/phages_plus_phage_genes.txt
	join -1 1 -2 1 <(cat /dev/shm/phages_that_contain_phage_genes.txt )  <(cat /dev/shm/all_phages_plus_phage_genes.txt) >> /dev/shm/phages_plus_phage_genes.txt
	## less /dev/shm/phages_plus_phage_genes.txt

	rm -f /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt
	sort -k2 /dev/shm/phages_plus_phage_genes.txt >> /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt

	##make a list of shared phage genes
	##the combinations are in phages_plus_phage_genes_sorted_on_genes.txt
	rm -f /dev/shm/shared_phage_genes.txt
	cat /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt | sort -k2 | uniq | tr " " "\t" | cut -f3 | uniq -c | sort -r -k 1 -n | 
	awk '{if($1!="1")print $2}' | sort >> /dev/shm/shared_phage_genes.txt
	## less /dev/shm/shared_phage_genes.txt

	##now we have the shorter list with gene content
	## join -1 1 -2 2 -o2.1,1.1 /dev/shm/shared_phage_genes.txt /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt| tr " " "\t" | less

	rm -f /dev/shm/phages_with_shared_pangenome.txt
	join -1 1 -2 2 -o2.1,1.1 /dev/shm/shared_phage_genes.txt /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt| tr " " "\t" | cut -f1| sort | uniq >> /dev/shm/phages_with_shared_pangenome.txt
	## less /dev/shm/phages_with_shared_pangenome.txt

	## and we are left with 566 (out of 583) phages that have some shared gene content from shared_phage_genes.txt
	##number of phages with shared pangenome
	echo "Number of phages with shared genes:"
	cat /dev/shm/phages_with_shared_pangenome.txt | wc -l

	##let us report the number of shared genes pangenome
	echo "Number of shared genes:"
	wc -l /dev/shm/shared_phage_genes.txt

	echo "Mostly shared genes (top $nr_of_genes_to_be_shown):"
	cat /dev/shm/phages_plus_phage_genes_sorted_on_genes.txt | sort -k2 | uniq | tr " " "\t" | cut -f3 | uniq -c | sort -r -k 1 -n \
		| awk '{if($1!="1")print $1,$2}' | head -${nr_of_genes_to_be_shown}


	##Now if we would sum $1 over the above lines and divide by number of phages and genes we have the coverage
}


function ana_iterations {

	##Output of VirSorter on 1000 samples with Viromes DB
	##out_dir=/hosts/linuxhome/mutant2/tmp/richard/virsorter_output/vs_all_1000
	out_dir=$1

	nr_of_genes_to_be_shown=20

	echo "Analyze VirSorter Output directory: $out_dir"

	nr_iterations=$(grep -c VIRSo -c $out_dir/r_*/global_signal_backup.csv | wc -l)

	nr_iterations=$((nr_iterations-1))

	echo "Nr. of iterations: $nr_iterations"

	i=0
	while [ $i -lt $nr_iterations ]
	do	i=$((i+1))
		iteration=r_$i
		analyze_one_iteration $out_dir $iteration $nr_of_genes_to_be_shown
	done
}
