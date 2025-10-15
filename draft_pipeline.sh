set -u

# Can do this at a higher level. Really just need beam parameters at sample, pdb file, and solvent information. 
# From there the below files can be generated. 

cd $(dirname "$0")

sim_handle=high_dmg
gromacs_unitcell_target_handle=4et8H_full_struct_Hfix  # TODO automate generation using "#%%--------STRUCTURE CONSTRUCTOR------" in scripts/scattering/scatter.py

TEMP_conversion_output_tag="_stoch" # TODO remove 


plasma_input_dir=input/pipeline_tests/
#targets_dir=scripts/scattering/targets/

gromacs_work_folder=~/cr-md/lys_pipeline_test/

mkdir -p $gromacs_work_folder


echo "Running plasma simulation"
./ac4dc_no_tbr $input_file_path

TEMP_serial_num=1 # TODO need to automate this.
sim_output_handle=${sim_handle}_${TEMP_serial_num}


echo "Preparing gromacs file"

bash pipeline_scripting/prepare_from_pdb.sh $gromacs_work_folder $gromacs_unitcell_target_handle



input_file_path=$plasma_input_dir/$sim_handle.mol
gromacs_file_path=$gromacs_work_folder/$gromacs_unitcell_target_handle.gro     # TODO automate generation

for expected_path in $input_file_path $gromacs_file_path; do 
    if [ ! -f $expected_path ]; then 
        echo "$expected_path not found" 
        exit 
    fi
done

echo "Generating molDStruct inputs from AC4DC"
python3.9 scripts/molDStructConversion/convert_to_molDStruct.py $sim_output_handle $gromacs_file_path
cp -r scripts/molDStructConversion/output/${sim_output_handle}${TEMP_conversion_output_tag}/IONIZATION_DATA $gromacs_work_folder



TEMP_nsteps=6001 # TODO read it
TEMP_dt=0.00002


echo Adding Lennard Jones params

bash pipeline_scripting/add_lennard_jones.sh $gromacs_work_folder 

echo "Running photon matter simulation"

bash pipeline_scripting/run_photon_matter_simulation.sh $gromacs_work_folder $gromacs_unitcell_target_handle.gro $TEMP_nsteps $TEMP_dt


MD_output_dir=scripts/scattering/targets/${sim_output_handle}_${gromacs_unitcell_target_handle}
mkdir -p MD_output_dir

bash pipeline_scripting/trjconv.sh $gromacs_work_folder $gromacs_file_path $MD_output_dir

python3.9 scripts/scattering/generate_MD_reflections.py $sim_output_handle $MD_output_dir

phenix.refine $OG_pdb_file_TODO $mtz_out_from_scatter_TODO





