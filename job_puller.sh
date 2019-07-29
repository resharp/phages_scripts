##examples

##get_job virsorter_scheduler.txt
##set_job_status virsorter_scheduler.txt "job_name" "status"

##job=$(get_job virsorter_scheduler.txt);echo $job


## get the first free job
function get_job {
	job_file=$1

	first_job=$(cat $job_file | grep -v "#" | grep -v "RUNNING"| grep -v "READY" | cut -f1 | head -1)

	echo $first_job
}

## get last job running
function last_job {
	job_file=$1

	last_job=$(cat $job_file | grep -v "#" | grep "RUNNING"| tail -1)

	echo $last_job
}

## return free_jobs
function view_free_jobs {

	job_file=$1

	cat $job_file | grep -v "#" | grep -v "RUNNING" | grep -v "READY" | cut -f1
}

function view_running {

	job_file=$1

	grep "RUNNING" $job_file | grep -v "READY" | cut -f1
}

function view_ready {
	job_file=$1

	grep "READY" $job_file | cut -f1
}

function set_job_status {
	job_file=$1
	job=$2
	status=$3

	##NB: double quotes matter here!
	##It appends the new status after the job, without erasing earlier statuses
	sed -i -r "s/($job)/\1\t$status/g" $job_file
}


function reset_jobs {

	job_file=$1
	rm -f /dev/shm/temp_job_file
	cat <(head -1 $job_file ) <(cat $job_file | grep -v "#" | cut -f1) >> /dev/shm/temp_job_file
	cat /dev/shm/temp_job_file > $job_file
}

function is_free_job_available {

	job_file=$1
	job=$(get_job $job_file)
	if [ $job ]; then echo "TRUE"; else echo "FALSE"; fi
}

