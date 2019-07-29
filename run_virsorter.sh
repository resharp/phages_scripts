##preparing the run script
source ~/scripts/job_puller.sh

##setting parameters (what mutant and the nr of samples)

##what is my hostname
mutant=$(hostname | tr '.' '\t' | cut -f1)

schedule_file=$1
job=$(get_job $schedule_file)

echo "---------------------------------------------------------------"
echo "Starting to run job run_virsorter on job $job on mutant $mutant"
echo "---------------------------------------------------------------"

##
set_job_status $schedule_file $job "RUNNING"
set_job_status $schedule_file $job $mutant

##take the samples from the job file
samples=$(cat ~/virsorter_admin/$job | cut -f1)

##copy the VirSorter database to the mutant if it does not exist
if [ ! -d "/hosts/linuxhome/$mutant/tmp/richard/virsorter-data" ]
then
	echo "starting to copy virsorter-data because they do not exist on $mutant yet"
	cp -r ~/virsorter-data/. /hosts/linuxhome/${mutant}/tmp/richard/virsorter-data
fi

echo "The job directory will be: "
job_dir=$(echo $job | tr '.' '_')
echo $job_dir

echo "Sample dir:" ##TODO mkdir
owner_dir=/hosts/linuxhome/${mutant}/tmp/richard
mkdir $owner_dir
sample_parent=/hosts/linuxhome/${mutant}/tmp/richard/patric_samples
mkdir $sample_parent
sample_dir=/hosts/linuxhome/${mutant}/tmp/richard/patric_samples/$job_dir
mkdir $sample_dir
echo $sample_dir

echo "Output will be put in: "
out_dir=/linuxhome/tmp/richard/virsorter_output/$job_dir
echo $out_dir

##copy the genomes to the mutant
for sample in $samples
	##do echo "copy ~/patric/patricdb_20190321/ftp.patricbrc.org/genomes/$sample/$sample.fna to $sample_dir"
	do cp /hosts/linuxhome/mgx/DB/PATRIC/patric/ftp.patricbrc.org/genomes/$sample/$sample.fna $sample_dir
done

##concatenate 
sample=all_1000
cat $sample_dir/*.fna >> $sample_dir/$sample.fna

echo "Running VirSorter on the concatenated fna file"
##run
## run on --db 1 ! on 8 cores
wrapper_phage_contigs_sorter_iPlant.pl -f $sample_dir/$sample.fna --db 1 --wdir $out_dir --ncpu 4 --data-dir /linuxhome/tmp/richard/virsorter-data

set_job_status $schedule_file $job "READY"
