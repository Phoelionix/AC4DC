import sys
import numpy as np
import argparse
import os.path as path
import os

#TODO need to adjust num time steps and grid update period based on fwhm

# parser = argparse.ArgumentParser(description = "Empty Description")
# parser.add_argument("infile_path")
# parser.add_argument("outfile_directory")
# parser.add_argument("-f", "--fast", help = "Uses shell approximation for atoms", required = False, default = False)
# parser.add_argument("-md", "--MD", help = "Generate MD file", required = False, default = False)
# parser.add_argument("-n", "--name", help = "Custom name for generated .mol file (no ext)", required = False, default = '')
# input_argument = parser.parse_args()

def main():
  # https://www.xfel.eu/sites/sites_custom/site_xfel/content/e35165/e46561/e46876/e179573/e179574/xfel_file179576/19042023_Parameters_factsheet_2024-01_Final_eng.pdf
  #ENERGIES = [6000,9000,12000]  # keV, 6-12keV is general operating range of EXFEL, 15keV is max.~ 9300 comes up in papers/specs repeatedly. 
  # FWHMS = [5,15,25]  # fs, Prefer to explore outside capability just for now, snce this will be biggest effect on sim time.
  # PHOTON_COUNTS = [0.01,0.1,1,10]   # 10^12 per square micrometre. 10^13 an approx. upper bound on EXFEL capabilities (though really 10^12 is probably closer to what is achievable).
 
  # ENERGIES = [12000]   
  # FWHMS = [5,10,15,25,50,100] 
  # PHOTON_COUNTS = [0.1,1,10,100,1000]  
  # #'''
  # ENERGIES = [6000,9000,12000,15000]
  # FWHMS = [5,10,25,50,100] 
  # PHOTON_COUNTS = [10]  
  # #'''
  # ENERGIES = [6000,9000,12000,15000]
  # FWHMS = [25] 
  # PHOTON_COUNTS = [0.1,1,10,100,1000]  
  ##
  #ENERGIES = [6000,9000,12000,15000,18000]
  #FWHMS = [7.5,10,15,25,50,100] 
  #PHOTON_COUNTS = [0.01,0.1,1,10,100]  
  ##

  ENERGIES = [6000,8000,10000,12000]
  FWHMS = [10,15,25] 
  PHOTON_COUNTS = [0.1,1,10]  


  batch_folder_parent_dir = "input/"

  if len(sys.argv) < 2:
    print('Usage: python3 generate_batch.py /path/to/file/to/copy.mol')
    return 1
  
  if sys.argv[1].partition('.')[-1] != 'mol':
    print('Invalid file: argument does not end in ".mol".')
    return 1

  # User provides path of input file
  infile_path = sys.argv[1] 
  handle = ''.join(path.basename(infile_path).partition('.')[:-1])[:-1]
  # We create the batch folder
  batch_dir = batch_folder_parent_dir + "batch_" + handle + "/"
  if not path.exists(batch_dir):
    os.makedirs(batch_dir)

  atoms, volume = copy_params(infile_path,handle)
  i = 1
  for fwhm in FWHMS:
    for photon_count in PHOTON_COUNTS:
      for energy in ENERGIES:
        feedin_dict = {
          'atoms' : atoms,
          'volume' : volume,
          'energy' : energy,
          'fwhm' : fwhm,
          'photon_count' : photon_count,
        }    
        outfile_path = None
        while outfile_path == None or path.exists(outfile_path):
          outfile_name = handle + "-" + str(i)
          outfile_path = batch_dir + outfile_name + '.mol'
          i+=1
        make_mol_file(infile_path, outfile_path,**feedin_dict)

  return 0

def copy_params(mol_path,handle):
    '''
    Reads the control file and returns the relevant parameters within
    '''    
    input_path = "input/"
    molfile = find_mol_file_from_directory(input_path,handle)
    atom_lines = []
    volume = None    
    with open(molfile, 'r') as f:
        reading_atoms = False
        reading_volume = False
        n = 0
        for line in f:
            if line.startswith("#ATOMS"):
                reading_atoms = True
                continue                
            elif line.startswith("#VOLUME"):
                reading_volume=True
                continue
            elif line.startswith("#") or line.startswith("//") or len(line.strip()) == 0:
                reading_atoms = False
                reading_volume = False
                n = 0
                continue
            if line.startswith("####END####"):
                break
            if reading_atoms:
                atom_lines.append(line)
            if reading_volume:
                if n == 0:
                    volume = float(line.split(' ')[0])
                n += 1
    print("Atoms:",atom_lines)
    print("Volume:", volume)
    return atom_lines, volume

def make_mol_file(fname, outfile, **incoming_params):
  print("Outputting to " + outfile)
  atoms = incoming_params.get('atoms')
  volume = incoming_params.get('volume')
  energy = incoming_params.get('energy')
  fwhm = incoming_params.get('fwhm')
  photon_count = incoming_params.get('photon_count')

  update_period = min(10,fwhm/3)
  num_steps = max(1000,round(fwhm*50))

  # Creare plasma input file
  plasma_file = open(outfile, 'w')
  plasma_file.write("""// AC4DC input file generated from: """ + fname + """\n\n""")
  plasma_file.write("""#ATOMS\n""")
  for line in atoms:
    plasma_file.write(line)
  
  plasma_file.write("""\n#VOLUME\n""")
  plasma_file.write("""%.2f      // Volume per molecule in Angstrom^3.\n""" % volume)
  plasma_file.write("""2000         // Radius of a sample in Angstrom. Used for effective escape rate of photo-electrons.\n""")
  plasma_file.write("""none         // Spatial boundary shape - options are none (confined system), spherical, cylindrical, planar.\n""")

  plasma_file.write("""\n#PULSE\n""")
  plasma_file.write("""%.0f         // Photon energy in eV.\n""" % energy)
  plasma_file.write("""%.1f         // Pulse-dependent definition of duration in femtoseconds: FWHM for Gaussian pulse, pulse width for square pulse.\n""" % fwhm)
  plasma_file.write("""gaussian     // Pulse shape.\n""")

  plasma_file.write("""\n#USE_COUNT\n""")
  plasma_file.write("""true         // enabled? true/false.\n""")
  plasma_file.write("""%.2f            // Photon density (x10^12 ph.µm-2)\n""" %photon_count)

  plasma_file.write("""\n#NUMERICAL\n""")
  plasma_file.write("""%.0f        // Initial guess for number of time step points for rate equation.\n""" %num_steps)
  plasma_file.write("""18           // Number of threads in OpenMP.\n""")

  plasma_file.write("""\n#DYNAMIC_GRID\n""")
  plasma_file.write("""medium       // Grid regions preset, options are 'low', 'medium', 'high', (accuracy) among others (see README).\n""")
  plasma_file.write("""%.2f         // Grid update period in fs, (dynamic grid only).\n""" %update_period)

  plasma_file.write("""\n#OUTPUT\n""")
  plasma_file.write("""800          // Number of time steps in the output files.\n""")
  plasma_file.write("""4000         // Number of free-electron grid points in output file.\n""")

  plasma_file.write("""\n#DEBUG\n""")
  plasma_file.write("""999           // Simulation cutoff time in fs (can end before but not after).\n""")
  plasma_file.write("""0.01          // Interval to display current timestep [fs].\n""")
  plasma_file.write("####END####\n\n")
  
  plasma_file.write("Notes:")
  
  plasma_file.close()


# copied from elsewhere. TODO give this its own file or something and import only in the function to stop HPC form being upset. 
def find_mol_file_from_directory(input_directory, mol):
        # Get molfile from all subdirectories in input folder.
        molfname_candidates = []
        for dirpath, dirnames, fnames in os.walk(input_directory):
            #if not "input/" in dirpath: continue
            for molfname in [f for f in fnames if f == mol+".mol"]:
                molfname_candidates.append(path.join(dirpath, molfname))
        # No file found
        if len(molfname_candidates) == 0:
            print('\033[91m[ Error: File not found ]\033[0m ' + "No '"+ mol + ".mol' input file found in input folders, trying default path anyway.")
            molfile = input_directory+mol+".mol"
        elif len(molfname_candidates) == 1:
            molfile = molfname_candidates[0]
        # File found in multiple directories
        else:
            print('\033[95m' + "Multiple mol files with given name detected. Please input number corresponding to desired directory." + '\033[0m' )
            for idx, val in enumerate(molfname_candidates):
                print(idx,val)
            selected_file = False
            # Loop until user confirms file
            while selected_file == False: 
                molfile = molfname_candidates[int(input("Input directory number: "))]
                print('\033[95m' + molfile + " selected." + '\033[0m')
                y_n_check = input("Input 'y'/'n' to continue/select a different file: ")
                if y_n_check.casefold() not in map(str.casefold,["y","yes"]):
                    if y_n_check.casefold() not in map(str.casefold,["n","no"]):
                        print("Unknown response, using 'n'.")
                    continue
                else:
                    print("Continuing...")
                    selected_file = True
        return molfile    

if __name__ == "__main__":
  main()