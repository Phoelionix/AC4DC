set -u

gmx='/home/speno/programs/bin' # path to where GROMACS is installed


working_folder=$1
NUM=0

cd $(dirname $0)

rm -f gmxdump_out

"$gmx/gmxdump" -s $working_folder/output_$NUM.tpr >> gmxdump_out
#"$gmx/gmxdump" -p /home/speno/cr-md/$working_folder/output.top >> out