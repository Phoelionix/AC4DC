set -u

gmx='/home/speno/programs/bin' # path to where GROMACS is installed

#http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html

# BROAD outline of pipeline, still quite manual as need to fill in file params in each step currently.
# Step 1 run AC4DC sim
# Step 2 construct structure with symmetries (e.g. see "#%%--------STRUCTURE CONSTRUCTOR------" in scripts/scattering/scatter.py)
# Step 3 run this file
# Step 4 run scripts/molDStructConversion/convert_to_molDStruct.py and add necessary files (also lennard-jones)
# step 4.5 YOU MUST RUN bash gmxdump.sh and then python3 create_lennard_jones_file.py and move the lennard jones file into IONIZATION_DATA!!!  
# ^ There will be slight order changes in the lj file, possibly due to atom ordering. Just make sure you do this or you will be very confused (see dreaded "Negative" error) followed by "md->massT[ai]=..."
# Step 5 bash run_photon_matter_simulation.sh

# working_folder=Lys_salt
# target_handle=4et8H_full_struct_Hfix

working_folder=$1 # full path
target_handle=$2
a=$3
b=$4
c=$5
add_disordered_water=$6

working_folder="$(realpath $working_folder)"

mkdir -p $working_folder/_bash
cd $working_folder/_bash


name_of_solvate_program=genbox # not solvate in older version...

rm -f -r topology

mkdir topology

out1=$working_folder/_bash/${target_handle}_pre_water_processed.gro
rm -f $out1
#$gmx/pdb2gmx -f /home/speno/AC4DC/scripts/scattering/targets/${target_handle}.pdb -o $out1  -p  $working_folder/_bash/topol.top  -water spce  <<EOF
$gmx/pdb2gmx -f /home/speno/AC4DC/scripts/scattering/targets/${target_handle}.pdb -o $out1  -p  $working_folder/_bash/topol.top  -water spce  -nocmap <<EOF
9
EOF
# 9  CHARMM36 (2020). (works for HISD, HIS1, HEME)
# 8 CHARMM27
# EOF
# 15 OPLS-AA/L
# EOF


for expected_path in $out1; do 
    if [ ! -f $expected_path ]; then 
        echo "Error: $expected_path not found" 
        exit 
    fi
done

# https://pmc.ncbi.nlm.nih.gov/articles/PMC11457149/

out2=$out1
if $add_disordered_water; then 
    out2=$working_folder/_bash/${target_handle}_processed.gro
    rm -f $out2
    rm -f box.gro
    $gmx/editconf -f $out1 -o box.gro -c  -box $a $b $c # 7.9 7.9 3.8

    #NOte the p flag is to update the topology
    $gmx/$name_of_solvate_program -cp box.gro -cs spc216.gro -o $out2 -p topol.top
fi

for expected_path in $out2; do 
    if [ ! -f $expected_path ]; then 
        echo "Error: $expected_path not found" 
        exit 
    fi
done


######$gmx/$name_of_solvate_program -cp $out1 -cs spc216.gro -o $out2 


cp $out2  /home/speno/AC4DC/scripts/scattering/targets/${target_handle}.gro  # For AC4DC to add charges
cp $out2 $working_folder/${target_handle}.gro # For MD sim

$gmx/make_ndx -f $working_folder/${target_handle}.gro -o $working_folder/index.ndx << EOF
q
EOF

mv *.top topology/
mv *.itp topology/

if [ -d $working_folder/topology/ ]; then 
    rm -r $working_folder/topology/
fi


 
cp -r topology/ $working_folder/topology/ 




# #  w/o water  #
# grep -v HOH /home/speno/AC4DC/scripts/scattering/targets/4et8_constructed_struct.pdb > /home/speno/cr-md/$working_folder/4et8_clean.pdb