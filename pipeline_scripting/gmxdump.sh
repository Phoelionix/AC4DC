set -u

gmx='/home/speno/programs/bin' # path to where GROMACS is installed


working_folder=$1
gro_file=$2

cd $(dirname $0)

rm -f gmxdump_out


"$gmx/grompp"  -f $working_folder/full_sim.mdp -c $gro_file  -n $working_folder/index.ndx -p $working_folder/topology/topol.top -o $working_folder/output.tpr -maxwarn 3  # create .tpr file

rm mdout.mdp

"$gmx/gmxdump" -s $working_folder/output.tpr >> gmxdump_out
#"$gmx/gmxdump" -p $working_folder/topology/topol.top >> gmxdump_out

