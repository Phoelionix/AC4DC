set -x
#!/usr/bin/env


#ideal_handle=("ideal_group")
ideal_handle=("idealB")

num_ideal_results_to_sample=1

scattering_targets_folder=/home/speno/AC4DC/scripts/scattering/targets
scattering_output_folder=/home/speno/AC4DC/scripts/scattering/random_sample
ccp4_folder=/home/speno/CCP4_Workspace/lys_random_sample
refinement_model_tmp=model_solutions/$refinement_tmp.pdb

#comparison_model=$scattering_targets_folder/4et8_error.pdb
comparison_model=$scattering_targets_folder/4et8_major_error.pdb



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

refmac5 xyzin $comparison_model xyzout $refinement_model_tmp hklin $ideal_mtz_file hklout refmac_out.mtz  << eof-refmac 
NCYCLES 0
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

rm $refinement_model_tmp


# sfcheck hklin $mtz_file xyzin $refinement_model
# sfcheck hklin $mtz_file xyzin $target_model
# sfcheck hklin $ideal_mtz_file xyzin $target_model
# sfcheck hklin $ideal_mtz_file xyzin $refinement_model

