##################################################################################################################
## run_parallel_vs.sh
##
## Richard Gremmen
## gremmen@resharp.nl
##
## dependecies:
## bin/run_vs_1.sh	this is the script that runs VirSorter on one genome
## bin/nodeslist	list of servers to distribute jobs over
##
## script inspired by:
## https://spectraldifferences.wordpress.com/2015/04/26/execute-commands-on-multiple-computers-using-gnu-parallel-setting-up-a-cluster-on-the-cheap/
## https://www.slashroot.in/how-run-multiple-commands-parallel-linux
## e.g.
## parallel -a $par_input -j 5 echo "Number {}: Running on $(hostname)"
##
##
##################################################################################################################

#this can be tested on the input file
#par_input=patric_samples/capsid_samples.txt
## TODO use input file

#head -6 patric_samples/capsid_samples.txt > patric_samples/capsid_samples_6.txt
#par_input=patric_samples/capsid_samples_6.txt

par_input=$1

# parallel --sshloginfile nodeslist echo "Number {}: Running on \`hostname\`" ::: 1 2 3 4

##here is an example of taking an input file and running on multiple cores
##there is no distribution over servers yet

#parallel -a $par_input -j 5 bash bin/run_vs_1.sh {}

#with distribution over servers
parallel -a $par_input -j 12 --sshloginfile bin/nodeslist bash bin/run_vs_1.sh {}
