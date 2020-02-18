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
