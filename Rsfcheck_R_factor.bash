#set -x
#!/usr/bin/env

#model_handle=ideal_016
#model_handle=salt_017

# model_handle=ideal_030
# pattern_handle=4et8_fmodel
# model_handle=4et8
# pattern_handle=4et8_fmodel

model_handle=ideal_073
#model_handle=dmg_074
pattern_handle=lys_undamaged

ccp4_folder=/home/speno/CCP4_Workspace/ComparePatterns

cd $ccp4_folder
mkdir -p R_logs
mkdir -p R_logs
rm R_logs/${model_handle}.log # Note removes old log


pattern=endPipelineData/${pattern_handle}.mtz
model=endPipelineData/${model_handle}.pdb
#model_pattern=endPipelineData/${model_handle}.mtz # from phenix..

refinement_model_tmp=tmp/R${model_handle}_tmp.pdb



sfcheck hklin $pattern xyzin $model