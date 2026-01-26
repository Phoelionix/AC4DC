set -u

gmx='/home/speno/programs/bin' # path to where GROMACS is installed

xtc_file=$1
ionization_data_dir=$2
out_folder=$3



#gro_file_handle=lys_example

i=1
out_subfolder="$out_folder/$i"
while [ -d $out_subfolder ]; do
    i=$((i+1))
    out_subfolder="$out_folder/$i"
    if [ "$i" -ge 9999 ]; then
        echo "Critical error: Couldn't find a filename that doesn't exist!"
        exit
    fi 
done

#mkdir 
mv $ionization_data_dir/ $out_subfolder/
mv $xtc_file $out_subfolder/output_$i.xtc