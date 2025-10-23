#set -x
#!/usr/bin/env

#model_handle=real_observed_019
#model_handle=salt_031 # .3029
#model_handle=ideal_030 # 0.2964
#model_handle=real_029 # 0.2733
#model_handle=salt_032 # 0.3075
#model_handle=ideal_033 # 0.3119
#model_handle=real_034 # 0.267
model_handle=ideal_sigma_diff_035 # 0.3020
#model_handle=salt_017

experiment_handle=4et8-sf 

ccp4_folder=/home/speno/CCP4_Workspace/ComparePatterns

cd $ccp4_folder
mkdir -p R_logs
mkdir -p R_logs
rm R_logs/${model_handle}.log # Note removes old log


experimental_pattern=endPipelineData/${experiment_handle}.mtz
model=endPipelineData/${model_handle}.pdb
model_pattern=endPipelineData/${model_handle}.mtz # from phenix..

refinement_model_tmp=tmp/R${model_handle}_tmp.pdb

refmac5 xyzin $model xyzout $refinement_model_tmp hklin $experimental_pattern hklout refmac_out.mtz  << eof-refmac 
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
#LABIN FP=Fobs SIGFP=Sigma-Fobs FREE=FreeR_flag
END
eof-refmac

rm $refinement_model_tmp

