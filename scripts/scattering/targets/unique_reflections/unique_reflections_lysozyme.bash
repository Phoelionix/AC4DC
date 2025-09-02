# Requires CCP4

set -u

cd $(dirname "$0")

high_res=1.5
name="lysozyme"

unique hklout x_unq.mtz <<eof-unique
TITLE  Unique data for $name
LABOUT  F=FP SIGF=SIGFP
SYMM P43212
RESOL ${high_res}
CELL 79.000   79.000   38.000
eof-unique

freerflag HKLIN x_unq.mtz HKLOUT x_unq2.mtz <<eof-freerflag
END
eof-freerflag


cad HKLIN1 x_unq2.mtz HKLIN2 x_unq.mtz HKLOUT unique_reflections.mtz << eof-cad
LABI FILE 1  E1=FreeR_flag
LABI FILE 2  ALLIN
END
eof-cad

mtz2various hklin unique_reflections.mtz hklout unique_reflections_${name}_${high_res}.hkl << eof-mtz2various
OUTPUT USER '(3I4,2F7.1,I4)'
MISS 0 
END
eof-mtz2various

#mtz2hkl -f unique_reflections.mtz

