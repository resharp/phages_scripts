##################################################################################################################
## run_vs_1.sh
##
## Richard Gremmen
## gremmen@resharp.nl
##
## this scripts runs virsorter on one genome
## it can be run in parallel (e.g. from run_parallel_vs.sh)
##
##
## dependencies of this script
##
## Virsorter is run within a conda environment
## it depends on the databases that have to be installed on each server
## conda activate virsorter (if run from parallel please activate in user's bash startup file)
##
## example:
## bash bin/run_vs_1.sh 1895771.3
##################################################################################################################

genome=$1

nr_cores=1

mutant=$(hostname | tr '.' '\t' | cut -f1)

##  automatically installs the dependency on the databases on the mutant that the job
function copy_prerequisite_virsorter_database {
	if [ ! -d "/hosts/linuxhome/$mutant/tmp/richard/virsorter-data" ]
	then
		echo "starting to copy virsorter-data because they do not exist on $mutant yet"
		cp -r ~/virsorter-data/. /hosts/linuxhome/${mutant}/tmp/richard/virsorter-data
	fi
}

##copy the VirSorter database to the mutant if it does not exist
## this only has to be done in the first run
##TODO: make optional (put in args) this has been commented out because it is already present on all servers that are used in the distribution
# copy_prerequisite_virsorter_database


date_start=$(date)
log_line="${date_start} Running genome ${genome} on ${mutant}"

echo $log_line

##TODO: change directory of logging
##command line parameter or default setting
echo $log_line >> bin/run_vs_1_log.txt


##dependency on genome
##that has to be copied rom the central mgx server with all the genomes
#echo "The job directory will be: "
job_dir=$genome
#echo $job_dir

#echo "Sample dir:" ##TODO mkdir
#TODO: reuse mkdir if not exist like in the code above
owner_dir=/hosts/linuxhome/${mutant}/tmp/richard
mkdir $owner_dir
sample_parent=/hosts/linuxhome/${mutant}/tmp/richard/patric_samples
mkdir $sample_parent
sample_dir=/hosts/linuxhome/${mutant}/tmp/richard/patric_samples/$job_dir
mkdir $sample_dir
#echo $sample_dir

#echo "Output will be put in: "
out_dir=/linuxhome/tmp/richard/virsorter_output/$job_dir
#echo $out_dir

virsorter_data=/linuxhome/tmp/richard/virsorter-data

##copy the genome to the mutant
cp /hosts/linuxhome/mgx/DB/PATRIC/patric/ftp.patricbrc.org/genomes/$genome/$genome.fna $sample_dir


##echo "Running VirSorter"
##run
# run on --db 1 ! on 1 core
wrapper_phage_contigs_sorter_iPlant.pl -f $sample_dir/$genome.fna --db 1 --wdir $out_dir --ncpu $nr_cores --data-dir $virsorter_data


date_end=$(date)
log_line="${date_end} Finished genome ${genome} on ${mutant}"

echo $log_line
echo $log_line >> bin/run_vs_1_log.txt
