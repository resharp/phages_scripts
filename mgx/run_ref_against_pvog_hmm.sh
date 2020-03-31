#
# /hosts/linuxhome/chaperone/tmp/richard/guerin_data/all_refs.proteins.faa
#
# annotation in tabular form:
#
# pvog_hmm_file=/hosts/linuxhome/mgx/nikos/project-001-virusfunction/data/pvogs/all/all.hmm
# ref_dir=/hosts/linuxhome/chaperone/tmp/richard/guerin_data

function search_all_genes_all_pvogs {

	conda deactivate
	conda activate mgx

	#now we want to search all ref genomes for all
	ref_dir=$1
	pvog_hmm_file=$2

	profile=$pvog_hmm_file

	protein_file=${ref_dir}/all_refs.proteins.faa

	out_file=${ref_dir}/all_refs_pvog.txt
	out_table=${ref_dir}/all_refs_pvog_table.txt

	# search all genes in every genome against all hmm profiles
	echo "profiles: " $profiles
	echo "protein_file:" $protein_file

	time hmmsearch --cpu 12 --tblout $out_table $profile $protein_file > $out_file

}

