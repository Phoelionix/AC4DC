import os
import os.path as path
import numpy as np

def get_sim_params(input_path,molecular_path,handle):
    '''
    Reads the control file and returns the relevant parameters within
    '''    
    molfile = get_mol_file(input_path,molecular_path,handle,"y")
    outDir = molecular_path + handle 
    intFile = outDir + "/intensity.csv"
    raw = np.genfromtxt(intFile, comments='#', dtype=np.float64)
    #intensityData = raw[:,1]
    timeData = raw[:, 0]    
    start_t = timeData[0]; end_t = timeData[-1]
    with open(molfile, 'r') as f:
        reading_count = False
        reading_pulse = False
        n = 0
        for line in f:
            if line.startswith("#PULSE"):
                reading_pulse=True
                continue                
            elif line.startswith("#USE_COUNT"):
                reading_count=True
                continue
            elif line.startswith("#") or line.startswith("//") or len(line.strip()) == 0:
                reading_count = False
                reading_pulse = False
                n = 0
                continue
            if line.startswith("####END####"):
                break
            if reading_pulse:
                if n < 2:
                    val = float(line.split(' ')[0])         
                if n == 0:
                    energy = val
                elif n==1:
                    fwhm = val
                n += 1
            if reading_count:
                if n == 1:
                    incoming_photon_count = float(line.split(' ')[0])
                n += 1
    print("Time range:",start_t,"-",end_t)
    print("Energy:", energy)
    param_dict = ["Energy","Width","Count"]
    unit_dict = [" eV"," fs","×10^12"]
    return start_t,end_t, energy, fwhm, incoming_photon_count, param_dict,unit_dict 

def get_mol_file(input_path, molecular_path, mol, output_mol_query = ""):
    #### Inputs ####
    #Check if .mol file in outputs
    use_input_mol_file = False
    output_folder  = molecular_path + mol + '/'
    if not os.path.isdir(output_folder):
        raise Exception("\033[91m Cannot find simulation output folder '" + mol + "'\033[0m" + "(In directory: " +output_folder + ")" )
    molfile = output_folder + mol + '.mol'
    # Use same-named mol file in output folder by default
    if path.isfile(molfile):
            print("Mol file found in output folder " + output_folder)
    # Use any mol file in output folder, (allowing for changing folder name).
    else: 
        molfile = ""
        for file in os.listdir(output_folder):
            if file.endswith(".mol"):
                if molfile != "": 
                    molfile = ""
                    print("\033[91m[ Warning: File Ambiguity ]\033[0m Multiple **.mol files in output folder.")
                    break
                molfile = os.path.join(output_folder, file)
        if molfile != "":
            print(".mol file found in output folder " + output_folder)
    
    if path.isfile(molfile):
        y_n_index = 2 
        if output_mol_query != "":
            y_n_check = output_mol_query
            unknown_response_qualifier = "argument \033[92m'"   + output_mol_query + "'\033[0m not understood, using default .mol file."
        else:
            print('\033[95m' + "Use \033[94m'" + os.path.basename(molfile) +"'\033[95m found in output? ('y' - yes, 'n' - use a mol file from input/ directory.)" + '\033[0m')
            y_n_check = input("Input: ")
            unknown_response_qualifier = "Response not understood" + ", using 'y'."
        if y_n_check.casefold() not in map(str.casefold,["n","no"]):  #  User didn't say to use input\
            if y_n_check.casefold() not in map(str.casefold,["y","yes"]):
                print(unknown_response_qualifier)
            if output_mol_query != "":
                print('\033[95m' + "Using \033[94m'" + os.path.basename(molfile) +"'\033[95m found in output." + '\033[0m')
        else:
            print("Using mol file from " + input_path)
            use_input_mol_file = True
    else: 
        print("\033[93m[ Missing Mol File ]\033[0m copy of mol file used to generate output not found, searching input/ directory.\033[0m'" )
        use_input_mol_file = True
    if use_input_mol_file:       
        molfile = find_mol_file_from_directory(input_path,mol) 
        print("Using: " + molfile)
    return molfile

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