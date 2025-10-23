set -x
#!/usr/bin/env


#handles=("ideal_group" "no_salt_group" "salt_group")
#handles=("idealB")
handles=("ideal_group")
#handles=("no_salt_group")
#handles=("salt_group")
#handles=("no_salt_group" )
#handles=( "salt_group")
#handles=( "idealB")

#ideal_handle=("ideal_group")
ideal_handle=("idealB")

num_results_to_sample=60
num_ideal_results_to_sample=1
num_loops=1

scattering_output_folder=/home/speno/AC4DC/scripts/scattering/random_sample
scattering_targets_folder=/home/speno/AC4DC/scripts/scattering/targets
ccp4_folder=/home/speno/CCP4_Workspace/lys_random_sample
#initial_model=$ccp4_folder/4et8_no_B.pdb
#initial_model=$scattering_targets_folder/4et8_major_error.pdb
initial_model=$ccp4_folder/4et8_incomplete.pdb
sequence_model=$ccp4_folder/4et8_no_B.pdb
target_model=$ccp4_folder/4et8_no_B.pdb

num_refmacs_cycles=10
num_pipeline_iters=1


# IDEAL ############
ideal_mtz_file=$ccp4_folder/${ideal_handle}.mtz

rm $ideal_mtz_file
rm $scattering_output_folder/scalepack/${handle}_0.sca 

python3.9  scripts/scattering/sample_scalepack.py $ideal_handle $num_ideal_results_to_sample "0";
for tag in $(eval echo {1..$num_loops}); 
do 
for handle in "${handles[@]}";
do
	mtz_file=$ccp4_folder/${handle}_$tag.mtz
	

	# Remove files so doesn't run if fail
	rm $scattering_output_folder/scalepack/${handle}_$tag.sca 
	rm $mtz_file 

	python3.9  scripts/scattering/sample_scalepack.py $handle $num_results_to_sample $tag;
    
done;
done;