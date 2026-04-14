set -u

# args: $1: path_to_starting model,  $2: folder of reflections to run refinement runs with    
#bash pipeline_scripting/Refinement/refinement_batch.sh scripts/scattering/targets/2QSP.pdb scalepack_MDS_All/

xyz_path=$(realpath "$1");          # Initial model
#hkl_path=$(realpath "$2"); shift 2  # Reflections
hkl_folder=$(realpath "$2"); shift 2  # Reflections


#resolution=None
resolution=1.5
#wc=1
wc=1

num_duplicate_runs=3

remove_old='false'

base_dir=$(realpath $(dirname "$0"))


i=0

tmp_refinement_dir=tmp_refinement_$resolution-norm_wc

#for hkl_path in $hkl_folder/*.sca; do
while ((i < num_duplicate_runs)); do 
    ((i++))
    for hkl_path in $hkl_folder/*.mtz; do
        cd $base_dir
        xyz_file=${xyz_path##*/}
        xyz_handle=${xyz_file%.*}
        hkl_file=${hkl_path##*/}
        hkl_handle=${hkl_file%.*}



        paramFile=nondefault_params.eff

        out_handle=${xyz_handle}-${hkl_handle}

        ####
        mkdir -p $tmp_refinement_dir

        if $remove_old; then
            rm -f $tmp_refinement_dir/$out_handle/$paramFile
            rm -rf  $tmp_refinement_dir/$out_handle/
        fi
        mkdir -p $tmp_refinement_dir/$out_handle/

        cp $paramFile $tmp_refinement_dir/$out_handle/
        cp $xyz_path $tmp_refinement_dir/$out_handle/${xyz_handle}.pdb
        cp $hkl_path $tmp_refinement_dir/$out_handle/${hkl_handle}.mtz

        ####
        cd $tmp_refinement_dir/$out_handle

        #hemoglobin_2QSP-1_All_real-5-5.pdb

        serial_idx=1
        while [[ -f  "${out_handle}_${serial_idx}.pdb" ||  -f "${serial_idx}.lock" ]]; do 
            serial_idx=$((serial_idx+1))
        done

        touch ${serial_idx}.lock
        
        phenix.refine ${xyz_handle}.pdb ${hkl_handle}.mtz $paramFile serial=$serial_idx xray_data.high_resolution=$resolution output.overwrite=True output.prefix=$out_handle xray_data.r_free_flags.generate=True wc=$wc

        rm ${serial_idx}.lock
        cd ../
    done
done