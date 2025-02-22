set -x
#!/usr/bin/env


handles=("ideal_group" "no_salt_group" "salt_group")
num_results_to_sample=30
num_loops=5

scattering_output_folder=/home/speno/AC4DC/scripts/scattering/random_sample
ccp4_folder=/home/speno/CCP4_Workspace/lys_random_sample
initial_model=$ccp4_folder/4et8_quite_incomplete.pdb
sequence_model=$ccp4_folder/4et8.pdb


for tag in $(eval echo {1..$num_loops}); 
do 
for handle in "${handles[@]}";
do
	mtz_file=$ccp4_folder/${handle}_$tag.mtz

	python3.9  scripts/scattering/sample_scalepack.py $handle $num_results_to_sample $tag;
	scalepack2mtz hklin $scattering_output_folder/scalepack/${handle}_$tag.sca hklout $mtz_file << eof-scalepack2mtz 
	END
eof-scalepack2mtz

	sftools << eof-sftools
read $mtz_file
sort h k l
I2F col IMEAN SIGIMEAN
write $mtz_file
Y
END
eof-sftools

	cd $ccp4_folder
	mkdir -p model_solutions
	mkdir -p buccaneer_logs
	rm buccaneer_logs/${handle}_$tag.log # Note removes old log

	buccaneer_pipeline -mtzin $mtz_file -pdbin $initial_model -seqin $sequence_model -colin-fo Fobs,Sigma-Fobs -pdbout model_solutions/${handle}_$tag.pdb >> buccaneer_logs/${handle}_$tag.log

done; 
done;

