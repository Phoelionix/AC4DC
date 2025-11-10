set -u
idx=$1
working_folder=$2 # full path
filename=$3
num_steps=$4
dt=$5
shake=$6

gmx='/home/speno/programs/bin' # path to where GROMACS is installed

shake_script=$(realpath $(dirname $0))/shake.py 

cd $(dirname $0)

cp template_run_file.mdp $working_folder/full_sim.mdp

cd $working_folder

sed "s/NSTEPS_PLACEHOLDER/${num_steps}/g" full_sim.mdp > tmp.$$
mv tmp.$$ full_sim.mdp

sed "s/DT_PLACEHOLDER/${dt}/g" full_sim.mdp > tmp.$$
mv tmp.$$ full_sim.mdp




shaken_file=4et8H_full_struct_Hfix_shaken_${idx}.gro 


#filename="4et8.gro"


python3.9 $shake_script $filename $shaken_file $shake


rm -f output_$idx.tpr
# Loop because fails sometimes, maybe due to LD random seed? Need to test.
i=0
while [ ! -f output_$idx.tpr ]; do
    "$gmx/grompp" -f full_sim.mdp -po mdout_$idx -c $shaken_file -n index.ndx -p ./topology/topol.top -o output_$idx.tpr -maxwarn 3  # create .tpr file
    i=$((i+1))
    if [ "$i" -ge 99 ]; then
        echo "Couldn't generate tpr file"
        exit
    fi 
done


rm -f \#*.*.*\#  # Remove extra backup files like "#output_8.edr.10#"

"$gmx/mdrun" -s output_$idx.tpr -deffnm output_$idx -v -nt 16 # run simulation, use -nt X, where X is number of cores you want to run with specific number of cores. -v is verbose

rm MPI_slice_n*



# REMOVES OLD OUTPUTS
#rm *#
