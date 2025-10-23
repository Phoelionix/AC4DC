set -x
#!/usr/bin/env


#handles=("ideal_group" "no_salt_group" "salt_group")
#handles=("idealB")
handles=("salt_group")
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
scalepack2mtz hklin $scattering_output_folder/scalepack/${ideal_handle}_0.sca hklout $ideal_mtz_file << eof-scalepack2mtz 
END
eof-scalepack2mtz

sftools << eof-sftools
read $ideal_mtz_file
sort h k l
I2F col IMEAN SIGIMEAN
write $ideal_mtz_file
Y
END
eof-sftools

freerflag hklin $ideal_mtz_file hklout ${ideal_mtz_file}_tmp << eof-freerflag
END
eof-freerflag
mv ${ideal_mtz_file}_tmp $ideal_mtz_file

####################

for tag in $(eval echo {1..$num_loops}); 
do 
for handle in "${handles[@]}";
do
	mtz_file=$ccp4_folder/${handle}_$tag.mtz
	

	# Remove files so doesn't run if fail
	rm $scattering_output_folder/scalepack/${handle}_$tag.sca 
	rm $mtz_file 

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


# END OF REFINEMENT 

# Compare final model stats

# refmac5 xyzin $refinement_model xyzout $refinement_model_tmp hklin $mtz_file hklout refmac_out.mtz  << eof-refmac 
# NCYCLES 0
# WEIGHT AUTO
# MAKE HYDR NO
# REFI BREF ISOT
# MAKE NEWLIGAND NOEXIT
# SCALE TYPE SIMPLE
# SOLVENT YES
# PHOUT
# MONI DIST 1000000
# PDBOUT KEEP USERS
# LABIN FP=Fobs SIGFP=Sigma-Fobs FREE=FreeR_flag
# END
# eof-refmac

# refmac5 xyzin $refinement_model xyzout $refinement_model_tmp hklin $ideal_mtz_file hklout refmac_out.mtz  << eof-refmac 
# NCYCLES 0
# WEIGHT AUTO
# MAKE HYDR NO
# REFI BREF ISOT
# MAKE NEWLIGAND NOEXIT
# SCALE TYPE SIMPLE
# SOLVENT YES
# PHOUT
# MONI DIST 1000000
# PDBOUT KEEP USERS
# LABIN FP=Fobs SIGFP=Sigma-Fobs FREE=FreeR_flag
# END
# eof-refmac


# refmac5 xyzin $target_model xyzout $refinement_model_tmp hklin $ideal_mtz_file hklout refmac_out.mtz  << eof-refmac 
# NCYCLES 0
# WEIGHT AUTO
# MAKE HYDR NO
# REFI BREF ISOT
# MAKE NEWLIGAND NOEXIT
# SCALE TYPE SIMPLE
# SOLVENT YES
# PHOUT
# MONI DIST 1000000
# PDBOUT KEEP USERS
# LABIN FP=Fobs SIGFP=Sigma-Fobs FREE=FreeR_flag
# END
# eof-refmac

# rm $refinement_model_tmp


# sfcheck hklin $mtz_file xyzin $refinement_model
# sfcheck hklin $mtz_file xyzin $target_model
# sfcheck hklin $ideal_mtz_file xyzin $target_model
# sfcheck hklin $ideal_mtz_file xyzin $refinement_model

done; 
done; 
done;

