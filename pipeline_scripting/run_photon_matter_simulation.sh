set -u
idx=$1
working_folder=$2 # full path
base_gro_file=$3
num_steps=$4
dt=$5
shake=$6
vacuum=$7

gmx='/home/speno/programs/bin' # path to where GROMACS is installed

shake_script=$(realpath $(dirname $0))/shake.py 

cd $(dirname $0)

### .mdp file ###
if [ -f "$working_folder/full_sim.mdp" ]; then
    mv $working_folder/full_sim.mdp $working_folder/full_sim.mdp#
fi
cp template_run_file.mdp $working_folder/full_sim.mdp

cd $working_folder

sed "s/NSTEPS_PLACEHOLDER/${num_steps}/g" full_sim.mdp > tmp.$$
mv tmp.$$ full_sim.mdp

sed "s/DT_PLACEHOLDER/${dt}/g" full_sim.mdp > tmp.$$
mv tmp.$$ full_sim.mdp

if $vacuum; then
    sed "s/pbc                      = xyz/pbc                      = no/g" full_sim.mdp > tmp.$$
    mv tmp.$$ full_sim.mdp
fi

### end .mdp file ###
# NOTE now in $working_folder


if (( shake > 0 )); then
    MD_input_gro_file=input_file_shaken${idx}.gro 
    python $shake_script $base_gro_file $MD_input_gro_file $shake
else
    MD_input_gro_file=input_file${idx}.gro
    cp $base_gro_file $MD_input_gro_file
fi


rm -f output_$idx.tpr
# Loop because fails sometimes, maybe due to LD random seed? Need to test.
i=0
while [ ! -f output_$idx.tpr ]; do
    rm -f mdout.mdp
    "$gmx/grompp" -f full_sim.mdp -po mdout -c $MD_input_gro_file -n index.ndx -p ./topology/topol.top -o output_$idx.tpr -maxwarn 3  # create .tpr file
    i=$((i+1))
    if [ "$i" -ge 2 ]; then
        echo "Couldn't generate tpr file"
        exit
    fi 
done


rm -f \#*.*.*\#  # Remove extra backup files like "#output_8.edr.10#"

#echo "Running:  "$gmx/mdrun" -s output_$idx.tpr -deffnm output_$idx -v -nt 24"
"$gmx/mdrun" -s output_$idx.tpr -deffnm output_$idx -v -nt 1  # run simulation, use -nt X, where X is number of cores you want to run with specific number of cores. -v is verbose

rm MPI_slice_n*



# REMOVES OLD OUTPUTS
#rm *#
