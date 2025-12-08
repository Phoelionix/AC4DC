set -ue

SCRIPT_DIR=$(
  CDPATH= cd -P -- "$(dirname -- "$BASH_SOURCE")" && pwd
)

cd $SCRIPT_DIR


num_traj_to_sample=5

gmx='/home/speno/programs/bin'




stem=hemoglobin_2QSP
#sample_interval=25 # 0.0004
sample_interval=200 # 25

#for sim_idx in 4 5 7 8; do 
for sim_idx in 9; do 

  ###
  #gromacs_unitcell_target_handle=hemoglobin_solv_Hfix
  #ordered_ground_truth_pdb=scripts/scattering/targets/hemoglobin_no_solv.pdb
  #ordered_ground_truth_pdb=$SCRIPT_DIR/scripts/scattering/targets/hemoglobin_no_solv_single_asu.pdb

  gromacs_unitcell_target_handle=2qsp_unit_Hfix
  ordered_ground_truth_pdb=$SCRIPT_DIR/scripts/scattering/targets/2qspH_zero_B.pdb


  gromacs_file_path=$SCRIPT_DIR/scripts/scattering/targets/$gromacs_unitcell_target_handle.gro
  #####

  # TODO use following to determine num frames and timestep to get fixed num frames?:
  #  ${gmx}check -f file.xtc &> log 


  plasma_handle=$stem-$sim_idx 
  MD_output_dir=scripts/scattering/targets/${plasma_handle}_${gromacs_unitcell_target_handle}/

  tmp_out_dir=$MD_output_dir/snapshots/

  rm -rf $tmp_out_dir
  mkdir $tmp_out_dir
  i=0
  for xtc_file in $MD_output_dir/*.xtc; do    
      bash pipeline_scripting/trjconv.sh $xtc_file $tmp_out_dir $gromacs_file_path $sample_interval
      i=$((i+1))
      if [ "$i" -ge $num_traj_to_sample ]; then
          break
      fi 
  done

  python3.9 scripts/scattering/generate_MD_reflections.py $plasma_handle $tmp_out_dir $ordered_ground_truth_pdb # $ordered_ground_truth_pdb

  #rm -r $tmp_out_dir

done