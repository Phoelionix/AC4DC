set -x -u
#!/usr/bin/env



folder=scripts/scattering/targets
xyzin=$folder/4et8H.pdb
xyzout=$folder/4et8H_zero_B.pdb



pdbset xyzin $xyzin xyzout $xyzout << eof 
BFACTOR ALWAYS 0
END
eof
