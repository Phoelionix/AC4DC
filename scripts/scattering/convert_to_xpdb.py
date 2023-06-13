# xpdb.py -- extensions to Bio.PDB
# (c) 2009 Oliver Beckstein
# Released under the same license as Biopython.
# See http://biopython.org/wiki/Reading_large_PDB_files
#%%
import sys
import Bio.PDB
import Bio.PDB.StructureBuilder
def main(infname, outfname):
  # Serial numbers are taken as remainders of 1e4
  folder = "solvate_1.0/"
  with open(folder+infname,'r') as infile, open(folder+outfname,'w') as outfile:
    line = infile.readline()
    while line!="":
        if line[0:4] == "ATOM" and line[0:6] != "ATOM  ":
            serial_number = line[4:].split()[0].strip()
            serial_number = str(int(int(serial_number)%1e5))
            line = "ATOM  "+ serial_number + line[6+len(serial_number):]
        outfile.write(line)
        line = infile.readline()


main("sol_4et8_8_cell.pdb","lys_8_cell.xpdb")
# %%
