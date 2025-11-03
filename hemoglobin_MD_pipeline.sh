set -u

# Can do this at a higher level. Really just need beam parameters at sample, pdb file, and solvent information. 
# From there the below files can be generated. 

NUM_MD=10

shake=0

TEMP_nsteps=2000 # TODO read it
TEMP_dt=0.00003

cd $(dirname "$0")

pdb_file="scripts/scattering/targets/hemoglobin_solv_Hfix.pdb" 

sim_handle=high_dmg
# NOTE had to change HSD to HIS to work with pdb2gmx ()
gromacs_unitcell_target_handle=hemoglobin_solv_Hfix  # TODO automate generation using "#%%--------STRUCTURE CONSTRUCTOR------" in scripts/scattering/scatter.py + phenix.ready_set for hydrogens

TEMP_conversion_output_tag="_stoch" # TODO remove 


plasma_input_dir=input/pipeline_tests/
#targets_dir=scripts/scattering/targets/

#gromacs_work_folder=~/cr-md/lys_pipeline_test/
gromacs_work_folder=~/cr-md/hemoglobin_uppsula/

mkdir -p $gromacs_work_folder


echo "Running plasma simulation"
#./ac4dc_no_tbr $input_file_path

TEMP_serial_num=1 # TODO need to automate this.
#sim_output_handle=${sim_handle}_${TEMP_serial_num}


echo "Preparing gromacs file"

#bash pipeline_scripting/prepare_from_pdb.sh $gromacs_work_folder $gromacs_unitcell_target_handle 6.251   8.097  11.148 false
# NOTE!!! manually removed the angles data from topol_Other_chain_W.itp after this... Bonds left unchanged... But why is this necessary? GROMACS version?


input_file_path=$plasma_input_dir/$sim_handle.mol
gromacs_file_path=$gromacs_work_folder/$gromacs_unitcell_target_handle.gro     # TODO automate generation

# for expected_path in $input_file_path $gromacs_file_path; do 
#     if [ ! -f $expected_path ]; then 
#         echo "$expected_path not found" 
#         exit 
#     fi
# done

    # echo "Generating molDStruct inputs from AC4DC"
    # #python3.9 scripts/molDStructConversion/convert_to_molDStruct.py $sim_output_handle $gromacs_file_path
    # cp -r scripts/molDStructConversion/output/${sim_output_handle}${TEMP_conversion_output_tag}/IONIZATION_DATA $gromacs_work_folder







new_charges_each_loop=false
if $new_charges_each_loop; then
    echo "Not implemented"
    exit
fi

echo "Running photon matter simulation"
width=10fs
for ((runid=1; runid<10; runid++)) {
    plasma_tag=${width}_${runid}
    moldstruct_conversion_output_handle=converted_charges_${plasma_tag}

    MD_output_dir=scripts/scattering/targets/${moldstruct_conversion_output_handle}_${gromacs_unitcell_target_handle}/
    #MD_output_dir=scripts/scattering/targets/
    mkdir -p $MD_output_dir
    
    python3.9 scripts/molDStructConversion/convert_h5.py scripts/molDStructConversion/hdf5_files/ $width $runid
    ION_DATA_PATH=scripts/molDStructConversion/output/${moldstruct_conversion_output_handle}${TEMP_conversion_output_tag}/IONIZATION_DATA
    
    if [ ! -d $ION_DATA_PATH ]; then 
        echo "$ION_DATA_PATH not found"
        exit
    fi
    
    for ((idx=1; idx<(NUM_MD+1); idx++)) {
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
        
        bash pipeline_scripting/trjconv.sh $idx $gromacs_work_folder $gromacs_file_path $MD_output_dir
    }
}





#python3.9 scripts/scattering/generate_MD_reflections.py $sim_output_handle $MD_output_dir

#phenix.refine $pdb_file $mtz_out_from_scatter_TODO





