set -u 

working_folder=$1


cd $(dirname $0)

bash gmxdump.sh $working_folder

python3.9 create_lennard_jones_file.py

out_lj_file=lennard_jones_parameters.txt


if [ ! -s $out_lj_file ]; then 
    echo Error: $out_lj_file is empty
    rm gmxdump_out  # gmxdump.sh dumps to 'out'
    exit 
fi



rm -f $working_folder/IONIZATION_DATA/$out_lj_file
mv $out_lj_file $working_folder/IONIZATION_DATA/


rm gmxdump_out
