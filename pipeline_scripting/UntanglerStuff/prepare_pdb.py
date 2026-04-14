from UntangleFunctions import prepare_pdb
import sys

prepare_pdb(*sys.argv[1:],allow_no_altloc=True)