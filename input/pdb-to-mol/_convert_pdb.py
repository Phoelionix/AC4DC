# WARNING seems like this isn't working.

import sys
import os
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "Empty Description")
parser.add_argument("infile_path")
parser.add_argument("outfile_directory")
parser.add_argument("-f", "--fast", help = "Uses shell approximation for atoms", required = False, default = False)
parser.add_argument("-md", "--MD", help = "Generate MD file", required = False, default = False)
parser.add_argument("-n", "--name", help = "Custom name for generated .mol file (no ext)", required = False, default = '')
input_argument = parser.parse_args()

def main():
  if len(sys.argv) < 3:
    print('Usage: converter path/to/pdb_files/protein.pdb path/to/mol_files/ path/to/outdir')
    return 1
  
  if sys.argv[1].partition('.')[-1] != 'pdb':
    print('No input .pdb file provided.')
    return 1

  # Get paths for in and out files, based on input arguments. 
  infile_path = sys.argv[1]
  outfile_dir = sys.argv[2] + ('/' if sys.argv[2][-1] != '/' else '') if sys.argv[2] !='/' else ''
  input_argument.name = '_' + input_argument.name 
  if input_argument.name != '_': input_argument.name += '_'
  outfile_name = 'pdb' + input_argument.name + os.path.splitext(infile_path)[0]
  outfile_path = outfile_dir + outfile_name + '.mol'
  out_MD_path = sys.argv[2] + 'MD_' + os.path.splitext(infile_path)[0] + '.xyz'

  # Read Harry Format pdb file
  atoms = read_pdb(infile_path)

  atom_cnt = 0
  # Output atomic composition.
  for k, v in atoms.items():
    atom_cnt += len(v)
    print("{} : {} atoms".format(k, len(v)))


  # Output unit cell dimensions.
  uc_bounds = []
  for i in range(3):
    tmp = [0, 0]
    for k in atoms:
      atoms[k] = np.array(atoms[k])
      tmp[0] = min(tmp[0], np.min(atoms[k][:, i]))
      tmp[1] = max(tmp[1], np.max(atoms[k][:, i]))

    uc_bounds.append(tmp)
    print("{}:  min {} ; max {}".format(i, uc_bounds[i][0], uc_bounds[i][1]))

  feedin_dict = {
    'boundary' : uc_bounds,
    'atoms' : atoms
  }

  make_plasma_input(os.path.splitext(infile_path)[0] + ".pdb", outfile_path,**feedin_dict)
  if input_argument.MD:
    make_MD_input(out_MD_path,**feedin_dict)

  return 0


def read_pdb(infile):
  atoms = {}
  uc_bound = [[-0.5*79.0, 0.5*79.0], [-0.5*79.0, 0.5*79.0], [-19.0, 19.0]]
  traslations = [v[1] - v[0] for v in uc_bound]
  atom_names = ["H", "C", "N", "O", "P", "S"]
  file = open(infile)

  in_header = True
  for i, line in enumerate(file):
    # Skip boilerplate at start
    if in_header:
      if  line[0:5] == "ATOM ":
        in_header = False
      else: continue

    tmp = line.split()
    """
    if tmp[0] == "TER":
      # Atoms already recorded are inside a unit cell.
      # Determine it's boundaries and drop atoms beyond it.
      for i in range(3):
        boundary = [0, 0]
        for k in atoms:
          axis = np.array(atoms[k])
          boundary[0] = min(boundary[0], np.min(axis[:, i]))
          boundary[1] = max(boundary[1], np.max(axis[:, i]))
        uc_bound.append(boundary)
    """

    if tmp[-1] not in atom_names: continue
    name = tmp[-1]
    x = float(tmp[-6])
    y = float(tmp[-5])
    z = float(tmp[-4])
    if len(uc_bound) > 0:
      if x < uc_bound[0][0] or x > uc_bound[0][1]: continue
      if y < uc_bound[1][0] or y > uc_bound[1][1]: continue
      if z < uc_bound[2][0] or z > uc_bound[2][1]: continue

    if name in atoms:
      atoms[name].append([x, y, z])
    else:
      atoms[name] = [[x, y, z]]

  file.close()
  return atoms

# (Arguably should be called make_mol_file)
def make_plasma_input(fname, outfile, **kwargs):
  extra_text = ''
  if input_argument.fast:
    extra_text += ', using atoms w/ shell approximations for orbitals'
  print("Outputting to " + outfile + extra_text)
  uc_bounds = kwargs.get('boundary')
  atoms = kwargs.get('atoms')

  # Creare plasma input file template from pdb.
  plasma_file = open(outfile, 'w')
  plasma_file.write("""// Template AC4DC input file created from """ + fname + """\n\n""")
  plasma_file.write("""#ATOMS\n""")
  for k, v in atoms.items():
    if input_argument.fast:
      k += '_fast'
    plasma_file.write("""{} {}\n""".format(k, len(v)))
  
  # Create a volume section.
  plasma_file.write("""\n#VOLUME\n""")
  tmp = 1
  for coord in uc_bounds: tmp *= coord[1] - coord[0]
  plasma_file.write("""%.2f      // Volume per molecule in Angstrom^3.\n""" % tmp)
  plasma_file.write("""2000      // Radius of a sample in Angstrom. Used for effective escape rate of photo-electrons.\n""")
  plasma_file.write("""none      // Spatial boundary shape - options are none (confined system), spherical, cylindrical, planar.\n""")


  # Fill in th rest with default values.
  plasma_file.write("""\n#PULSE\n""")
  plasma_file.write("""6000         // Photon energy in eV.\n""")
  plasma_file.write("""10           // Pulse-dependent definition of duration in femtoseconds: FWHM for Gaussian pulse, pulse width for square pulse.\n""")
  plasma_file.write("""2500         // Pulse fluence in 10^4 * J/cm^2.\n""")
  plasma_file.write("""square       // Pulse shape.\n""")

  plasma_file.write("""\n#NUMERICAL\n""")
  plasma_file.write("""400000       // Initial guess for number of time step points for rate equation.""")
  plasma_file.write(""" (TODO following is not working currently?: -S.P.) If the value is 0, program skips rate equation solving step.\n""")
  plasma_file.write("""12           // Number of threads in OpenMP.\n""")
  plasma_file.write("""0            // [only for nonthermal] minimum free-electron energy in eV.\n""")
  plasma_file.write("""10000        // [only for nonthermal] maximum free-electron energy in eV. Should not be near photon energy.\n""")
  plasma_file.write("""100          // total number of free-electron grid points.\n""")
  plasma_file.write("""powerlaw     // electron grid type.\n""")
  plasma_file.write("""35           // Number of "low-energy" grid points (hybrid & powerlaw grids only).\n""")
  plasma_file.write("""2000         // transition energy in eV.\n""")

  plasma_file.write("""\n#OUTPUT\n""")
  plasma_file.write("""800          // Number of time steps in the output files.\n""")
  plasma_file.write("""4000         // Number of free-electron grid points in output file.\n""")

  plasma_file.write("""\n#DEBUG\n""")
  plasma_file.write("""1         // Proportion of time steps to iterate through before stopping early.\n""")
  plasma_file.write("""0.001         // Interval to update current timestep [fs].\n""")
  plasma_file.write("####END####\n\n")
  
  plasma_file.write("Notes:")

  plasma_file.close()


def make_MD_input(fname,**kwargs):
  uc_bounds = kwargs.get('boundary')
  atoms = kwargs.get('atoms')

 # Creare plasma input file template from Harry's pdb.
  md_file = open(fname, 'w')
  num_atoms = 0
  for k, v in atoms.items():
    num_atoms += len(v)
  md_file.write("""{}\n""".format(num_atoms))
  for elem in uc_bounds:
    md_file.write("""%.3f   """ % (elem[1] - elem[0]))

  # Atoms.
  x_shift = uc_bounds[0][0]
  y_shift = uc_bounds[1][0]
  z_shift = uc_bounds[2][0]

  for name in atoms:
    for atom in atoms[name]:
      md_file.write("""\n{}  {:.3f}  {:.3f}  {:.3f}""".format(name, atom[0] - x_shift, atom[1] - y_shift, atom[2] - z_shift))



if __name__ == "__main__":
  main()