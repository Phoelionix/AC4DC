set -u
working_folder=$1 # full path
filename=$2
num_steps=$3
dt=$4


gmx='/home/speno/programs/bin' # path to where GROMACS is installed

shake_script=$(realpath $(dirname $0))/shake.py 

cd $(dirname $0)

cp template_run_file.mdp $working_folder/full_sim.mdp

cd $working_folder

sed "s/NSTEPS_PLACEHOLDER/${num_steps}/g" full_sim.mdp > tmp.$$
mv tmp.$$ full_sim.mdp

sed "s/DT_PLACEHOLDER/${dt}/g" full_sim.mdp > tmp.$$
mv tmp.$$ full_sim.mdp



for ((idx=0; idx<1; idx++)) {
    shaken_file=4et8H_full_struct_Hfix_shaken_${idx}.gro 


    #filename="4et8.gro"


    python3.9 $shake_script $filename $shaken_file 0.01

   
    $gmx/grompp -f full_sim.mdp -po mdout_$idx -c $shaken_file -n index.ndx -p ./topology/topol.top -o output_$idx.tpr -maxwarn 3  # create .tpr file

    $gmx/mdrun -s output_$idx.tpr -v -deffnm output_$idx # run simulation, use -nt X, where X is number of cores you want to run with specific number of cores

    rm MPI_slice_n*
}


# REMOVES OLD OUTPUTS
#rm *#
