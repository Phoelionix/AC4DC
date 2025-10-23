#set -x
#!/usr/bin/env

#model_handle=ideal_016
#model_handle=salt_017

model_handle=salt_032 # 
#   R-factor           : 0.3043
#   Correlation factor :   0.8586
#model_handle=ideal_033 # 
#   R-factor           : 0.3011
#   Correlation factor :   0.8605
#model_handle=real_034 # 
#   R-factor           : 0.2634
#   Correlation factor :   0.8934

experiment_handle=4et8-sf # validation file

ccp4_folder=/home/speno/CCP4_Workspace/ComparePatterns

cd $ccp4_folder
mkdir -p R_logs
mkdir -p R_logs
rm R_logs/${model_handle}.log # Note removes old log


experimental_pattern=endPipelineData/${experiment_handle}.mtz
model=endPipelineData/${model_handle}.pdb
model_pattern=endPipelineData/${model_handle}.mtz # from phenix..

refinement_model_tmp=tmp/R${model_handle}_tmp.pdb



sfcheck hklin $experimental_pattern xyzin $model