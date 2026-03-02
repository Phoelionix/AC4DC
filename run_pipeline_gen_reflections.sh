set -ue

SCRIPT_DIR=$(
  CDPATH= cd -P -- "$(dirname -- "$BASH_SOURCE")" && pwd
)

cd $SCRIPT_DIR


num_traj_to_sample=10
num_snapshots=10 # per trajectory 



dt_of_MD_sim=0.01 # TODO read this or num_frames from trj file

stem=hemoglobin_2QSP

for sim_idx in 1; do 

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
  MD_output_parent_dir=$SCRIPT_DIR/MD_output/${plasma_handle}_${gromacs_unitcell_target_handle}/

  read -r timespan < <(python3.9 scripts/print_sim_params.py $plasma_handle timespan)
  num_frames=`echo $timespan $dt_of_MD_sim  | awk '{print $1/$2}'` # duration/dt
  sample_interval=`echo $num_frames $num_snapshots | awk '{print $1/$2}'`   # num_frames/num_snapshots
  sample_interval=$( printf "%.0f" $sample_interval)


  # Trajectory to pdb file of snapshots
  i=0
  for old_snapshots in $MD_output_parent_dir/*/snapshots.pdb; do
    if [ -f $old_snapshots ]; then
     rm $old_snapshots
    fi
  done
  for subdir in $MD_output_parent_dir/*/; do
      echo $subdir
      xtc_file=false
      for tmp in $subdir/*.xtc; do
        if $xtc_file; then 
          echo "Error: More than 1 xtc file in $subdir"
          exit
        fi
        xtc_file=$tmp
      done
      if [ ! -f $xtc_file ]; then 
        continue
      fi 
      if [ ! -s $xtc_file ]; then 
        continue
      fi 
      nice -n 5 bash pipeline_scripting/trjconv.sh $xtc_file $subdir/snapshots.pdb $gromacs_file_path $sample_interval
      i=$((i+1))
      if [ "$i" -ge $num_traj_to_sample ]; then
          break
      fi 
  done

  
nice -n 5 python3.9 scripts/scattering/generate_MD_reflections.py $plasma_handle $MD_output_parent_dir $ordered_ground_truth_pdb # $ordered_ground_truth_pdb

done
wait