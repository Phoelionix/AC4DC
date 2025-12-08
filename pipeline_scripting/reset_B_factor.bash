set -x -u
#!/usr/bin/env



folder=~/AC4DC/scripts/scattering/targets
#pdb_handle=4et8H
pdb_handle=2qspH

xyzin=$folder/$pdb_handle.pdb
xyzout=$folder/${pdb_handle}_zero_B.pdb



pdbset xyzin $xyzin xyzout $xyzout << eof 
BFACTOR ALWAYS 0
END
eof
