set -u

# Can do this at a higher level. Really just need beam parameters at sample, pdb file, and solvent information. 
# From there the below files can be generated. 


taskid=$1
NUM_MD=$2



stem=hemoglobin_2QSP
num_sims=9 
sim_idx=$(( (($taskid-1) % $num_sims)+1 ))

shake=0

plasma_handle=$stem-$sim_idx 




copy_topology_from_template=true

new_charges_each_loop=true



cd $(dirname "$0")






# NOTE had to change HSD to HIS to work with pdb2gmx ()
gromacs_unitcell_target_handle=2qsp_unit_Hfix # TODO automate generation. E.g. using "#%%--------STRUCTURE CONSTRUCTOR------" in scripts/scattering/scatter.py + phenix.ready_set for hydrogens


plasma_input_dir=input/pipeline_tests/


gromacs_work_folder=~/cr-md/${plasma_handle}-t${taskid}/
mkdir -p $gromacs_work_folder



#sim_handle=high_dmg
#input_file_path=$plasma_input_dir/$sim_handle.mol
#echo "Running plasma simulation"
#./ac4dc $input_file_path

#template_folder=hemoglobin_uppsala_template
template_folder_name=2qsp_template

template_folder=~/cr-md/$template_folder_name
add_disordered_water=true


CELL_LENGTHS_NM="6.5033   7.8273   10.9085"

######
if ! $copy_topology_from_template; then 
    echo "Preparing MD structural input files from pdb"
    bash pipeline_scripting/prepare_from_pdb.sh $template_folder $gromacs_unitcell_target_handle $CELL_LENGTHS_NM $add_disordered_water # TODO READ THIS
    bash pipeline_scripting/delete_dihedrals.sh $template_folder
    # NOTE!!! manually removed the angles data from topol_Other_chain_W.itp after this... Bonds left unchanged... But why is this necessary? GROMACS version?
fi
    echo "Taking MD structural input files from template"
    rm -rf $gromacs_work_folder/topology/
    cp -r $template_folder/topology/ $gromacs_work_folder  # Could make the scripts read from topology in a different directory from $gromacs_work_folder, but this is easiest.


    if [ ! -d $gromacs_work_folder/topology/ ]; then
        echo $gromacs_work_folder/topology/ does not exist
        exit
    fi

    # NOTE full_sim.mdp is transferred purely for generating LJ parameters. It will be overwritten.
    for required_file in $gromacs_unitcell_target_handle.gro index.ndx; do
        if [ ! -f "$gromacs_work_folder/$required_file" ]; then
            cp $template_folder/$required_file $gromacs_work_folder
        fi
        if [ ! -f "$gromacs_work_folder/$required_file" ]; then
            echo "$gromacs_work_folder/$required_file" does not exist
            exit
        fi
    done
######






gromacs_file_path=$gromacs_work_folder/$gromacs_unitcell_target_handle.gro     # TODO automate generation



#read -r dt nsteps < <(python3.9 scripts/print_sim_params.py $plasma_handle dt_ps nsteps)

echo "Running MD simulations"
DEBUG=false
for ((idx=1; idx<(NUM_MD+1); idx++)) {
    moldstruct_conversion_output_handle=$plasma_handle

    MD_output_dir=scripts/scattering/targets/${moldstruct_conversion_output_handle}_${gromacs_unitcell_target_handle}/
    mkdir -p $MD_output_dir
    
    ION_DATA_PATH=scripts/molDStructConversion/output/${moldstruct_conversion_output_handle}/IONIZATION_DATA
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

        python_command="python3.9 scripts/molDStructConversion/convert_to_molDStruct.py \
            $plasma_handle $gromacs_file_path $DEBUG"
        #nsteps=$($python_command | awk '/Creating charge file with/{print $5}') 
        read -r nsteps dt < <($python_command | awk '/Creating charge file with/{print $5, $7}')
        echo "$nsteps steps, dt = $dt ps"
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
        
        dummy_mdp=$gromacs_work_folder/full_sim.mdp
        cp pipeline_scripting/dummy_lj.mdp $dummy_mdp
        bash pipeline_scripting/add_lennard_jones.sh $gromacs_work_folder $gromacs_file_path
        out_lj_file=$gromacs_work_folder/IONIZATION_DATA/lennard_jones_parameters.txt
        if [ ! -s $out_lj_file ]; then 
            echo Error in generating LJ file
            exit 
        fi
        rm $dummy_mdp

        if [ -f $gromacs_work_folder/output_$idx.xtc ]; then
            mv $gromacs_work_folder/output_$idx.xtc $gromacs_work_folder/output_$idx.xtc#
        fi

        bash pipeline_scripting/run_photon_matter_simulation.sh $idx $gromacs_work_folder $gromacs_unitcell_target_handle.gro $nsteps $dt $shake
        if [ ! -f $gromacs_work_folder/output_$idx.xtc ]; then
            echo "output failed"
            exit
        fi

        bash pipeline_scripting/trjcopy.sh $gromacs_work_folder/output_${idx}.xtc  $MD_output_dir
        #bash pipeline_scripting/trjconv.sh $gromacs_work_folder/output_${idx}.xtc $MD_output_dir $gromacs_file_path
        
        if $DEBUG; then
            echo "DEBUG: EXITING EARLY"; exit
        fi
    ####
}





#python3.9 scripts/scattering/generate_MD_reflections.py $sim_output_handle $MD_output_dir

#phenix.refine $pdb_file $mtz_out_from_scatter_TODO





