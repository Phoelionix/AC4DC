set -x
#!/usr/bin/env


handles=("ideal_group" "no_salt_group" "salt_group")
num_results_to_sample=60
num_loops=1

scattering_output_folder=/home/speno/AC4DC/scripts/scattering/random_sample
ccp4_folder=/home/speno/CCP4_Workspace/lys_random_sample
initial_model=$ccp4_folder/4et8.pdb
sequence_model=$ccp4_folder/4et8.pdb

num_refmacs_cycles=10
num_pipeline_iters=1


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

freerflag hklin $mtz_file hklout ${mtz_file}_tmp << eof-freerflag
END
eof-freerflag

mv ${mtz_file}_tmp $mtz_file

### Refinement pipeline, mimicing that of ccp4i2
cd $ccp4_folder
mkdir -p model_solutions
mkdir -p buccaneer_logs
rm buccaneer_logs/${handle}_$tag.log # Note removes old log

refinement_model=model_solutions/${handle}_$tag.pdb
refinement_model_tmp=model_solutions/${handle}_${tag}_tmp.pdb

csheetbend -pdbin $initial_model -mtzin $mtz_file -pdbout $refinement_model -colin-fo Fobs,Sigma-Fobs -colin-free FreeR_flag -cycles 12 -resolution-by-cycle 6.0,3.0 -coord -radius-scale 4.0 

for i in $(eval echo {1..$num_pipeline_iters}); do

refmac5 xyzin $refinement_model xyzout $refinement_model_tmp hklin $mtz_file hklout refmac_out.mtz  << eof-refmac 
NCYCLES 10
WEIGHT AUTO
MAKE HYDR NO
REFI BREF ISOT
MAKE NEWLIGAND NOEXIT
SCALE TYPE SIMPLE
SOLVENT YES
PHOUT
MONI DIST 1000000
PDBOUT KEEP USERS
LABIN FP=Fobs SIGFP=Sigma-Fobs FREE=FreeR_flag
END
eof-refmac

mv $refinement_model_tmp $refinement_model

done; 
done; 
done;

