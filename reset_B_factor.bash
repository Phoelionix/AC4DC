set -x
#!/usr/bin/env



ccp4_folder=/home/speno/CCP4_Workspace/lys_random_sample
xyzin=$ccp4_folder/4et8.pdb
xyzout=$ccp4_folder/4et8_no_B.pdb



pdbset xyzin $xyzin xyzout $xyzout << eof 
BFACTOR ALWAYS 0
END
eof
