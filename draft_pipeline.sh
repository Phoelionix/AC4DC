set -u

# Can do this at a higher level. Really just need beam parameters at sample, pdb file, and solvent information. 
# From there the below files can be generated. 


taskid=single_node
NUM_MD=30
runid_step_size=1
id_start=7
id_end=16


shake=0
TEMP_nsteps=2000 # TODO read it
fluence_idx=1

copy_topology_from_template=true

new_charges_each_loop=true



cd $(dirname "$0")

pdb_file="scripts/scattering/targets/hemoglobin_solv_Hfix.pdb" 

# NOTE had to change HSD to HIS to work with pdb2gmx ()
gromacs_unitcell_target_handle=hemoglobin_solv_Hfix  # TODO automate generation using "#%%--------STRUCTURE CONSTRUCTOR------" in scripts/scattering/scatter.py + phenix.ready_set for hydrogens

TEMP_conversion_output_tag="_stoch" # TODO remove 


plasma_input_dir=input/pipeline_tests/
#targets_dir=scripts/scattering/targets/

#gromacs_work_folder=~/cr-md/lys_pipeline_test/
gromacs_work_folder=~/cr-md/hemoglobin_uppsala-${taskid}/
mkdir -p $gromacs_work_folder



echo "Preparing input files from pdb"

if $copy_topology_from_template; then 
    rm -rf $gromacs_work_folder/topology/
    cp -r ~/cr-md/hemoglobin_uppsala_template/topology/ $gromacs_work_folder  # Could make the scripts read from topology in a different directory from $gromacs_work_folder, but this is easiest.


    if [ ! -d $gromacs_work_folder/topology/ ]; then
        echo $gromacs_work_folder/topology/ does not exist
        exit
    fi

    for required_file in $gromacs_unitcell_target_handle.gro index.ndx; do
        if [ ! -f "$gromacs_work_folder/$required_file" ]; then
            cp ~/cr-md/hemoglobin_uppsala_template/$required_file $gromacs_work_folder
        fi
        if [ ! -f "$gromacs_work_folder/$required_file" ]; then
            echo "$gromacs_work_folder/$required_file" does not exist
            exit
        fi
    done
else
    bash pipeline_scripting/prepare_from_pdb.sh $gromacs_work_folder $gromacs_unitcell_target_handle 6.251   8.097  11.148 false
    # NOTE!!! manually removed the angles data from topol_Other_chain_W.itp after this... Bonds left unchanged... But why is this necessary? GROMACS version?
fi



#sim_handle=high_dmg
#input_file_path=$plasma_input_dir/$sim_handle.mol
#echo "Running plasma simulation"
#./ac4dc $input_file_path



gromacs_file_path=$gromacs_work_folder/$gromacs_unitcell_target_handle.gro     # TODO automate generation





echo "Running MD simulation"
for ((idx=1; idx<(NUM_MD+1); idx++)) {
for ((fluence_idx=0; fluence_idx<2; fluence_idx+=1)) {
    if (( fluence_idx == 0 )); then
        TEMP_dt=0.00003; width=10fs
    elif (( fluence_idx == 1 )); then
        TEMP_dt=0.00001; width=3fs
    else 
        echo "Fluence idx of $fluence_idx is invalid"
    exit
    fi

    for ((runid=id_start; runid<=id_end; runid+=runid_step_size)) {
        plasma_tag=${width}_${runid}
        moldstruct_conversion_output_handle=converted_charges_${plasma_tag}

        MD_output_dir=scripts/scattering/targets/${moldstruct_conversion_output_handle}_${gromacs_unitcell_target_handle}/
        mkdir -p $MD_output_dir
        
        ION_DATA_PATH=scripts/molDStructConversion/output/${moldstruct_conversion_output_handle}${TEMP_conversion_output_tag}/IONIZATION_DATA
        if $new_charges_each_loop; then
            if [ -d $ION_DATA_PATH ]; then 
                if [ -d $ION_DATA_PATH# ]; then
                    rm -r  $ION_DATA_PATH#
                fi
                mv $ION_DATA_PATH $ION_DATA_PATH#
            fi
        fi
        if [ ! -d $ION_DATA_PATH ]; then 
            echo "Generating electron data"
            python scripts/molDStructConversion/convert_h5.py scripts/molDStructConversion/hdf5_files/ $width $runid
        fi
        if [ ! -d $ION_DATA_PATH ]; then 
            echo "$ION_DATA_PATH was not generated for unknown reason."
            exit
        fi
        ####
            echo "Generating molDStruct inputs from AC4DC"
            #
            # 
            if [ -d ${gromacs_work_folder}/IONIZATION_DATA/ ]; then 
                rm -rf ${gromacs_work_folder}/IONIZATION_DATA#/
                mv ${gromacs_work_folder}/IONIZATION_DATA/ ${gromacs_work_folder}/IONIZATION_DATA#/
            fi 
            cp -r $ION_DATA_PATH $gromacs_work_folder
            echo Adding Lennard Jones params
            bash pipeline_scripting/add_lennard_jones.sh $gromacs_work_folder $gromacs_file_path 

            if [ -f $gromacs_work_folder/output_$idx.xtc ]; then
                mv $gromacs_work_folder/output_$idx.xtc $gromacs_work_folder/output_$idx.xtc#
            fi

            bash pipeline_scripting/run_photon_matter_simulation.sh $idx $gromacs_work_folder $gromacs_unitcell_target_handle.gro $TEMP_nsteps $TEMP_dt $shake
            if [ ! -f $gromacs_work_folder/output_$idx.xtc ]; then
                echo "output failed"
                exit
            fi

            bash pipeline_scripting/trjcopy.sh $gromacs_work_folder/output_${idx}.xtc  $MD_output_dir
        ####
    }
}
}





#python scripts/scattering/generate_MD_reflections.py $sim_output_handle $MD_output_dir

#phenix.refine $pdb_file $mtz_out_from_scatter_TODO





