import os
import os.path as path
import numpy as np
import sys
import re
import csv
from QoL import set_highlighted_excepthook

def get_sim_params(handle,input_path=None,molecular_path=None):
    '''
    Reads the control file and returns the relevant parameters within
    By default use input_path = "input/"
    '''    
    if input_path is None:
        input_path = path.abspath(path.join(__file__ ,"../../input/")) + "/"
    if molecular_path is None:
        molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
        
    molfile = get_mol_file(input_path,molecular_path,handle,"y",out_prefix_text = "Reading simulation parameters from")
    outDir = molecular_path + handle 
    intFile = outDir + "/intensity.csv"
    raw = np.genfromtxt(intFile, comments='#', dtype=np.float64)
    #intensityData = raw[:,1]
    timeData = raw[:, 0]    
    start_t = timeData[0]; end_t = timeData[-1]
    # for convenience the input file has multiple options for determining fluence: "Count", "Fluence", or "Intensity"    
    active_photon_measure_found = False
    photon_measure = None  
    reading_photons = False
    reading_pulse = False    
    reading_electron_source = False
    photon_unit = None
    photon_measure_val = None
    source_energy = None
    source_fraction = None
    source_duration = None

    def init_photon_read(param_type,unit):
        nonlocal reading_photons,photon_measure,photon_unit
        if active_photon_measure_found != False:
            raise Exception("two photon measures detected")                
        reading_photons = True   
        photon_measure = param_type  
        photon_unit = unit   
    with open(molfile, 'r') as f:
        n = 0
        for line in f:
            if line.startswith("#PULSE"):
                reading_pulse=True
                continue                
            elif line.startswith("#USE_COUNT"): # (incoming photon count)
                init_photon_read("Count","×10^12")              
                continue
            elif line.startswith("#USE_FLUENCE"):
                init_photon_read("Fluence","×10^4 J/cm^2")              
                continue
            elif line.startswith("#USE_INTENSITY"):
                init_photon_read("Intensity","×10^19 W/cm^2")              
                continue
            elif line.startswith("#ELECTRON_SOURCE"):
                reading_electron_source = True
                continue                        
            elif line.startswith("#") or line.startswith("//") or len(line.strip()) == 0:
                reading_photons = False
                reading_pulse = False
                reading_electron_source = False 
                n = 0
                continue
            if line.startswith("####END####"):
                break
            if reading_pulse:
                if n < 2:
                    val = float(line.split(' ')[0])         
                if n == 0:
                    photon_energy = val
                elif n==1:
                    fwhm = val
                n += 1
            if reading_photons:
                if n == 0:
                    if (line.split(' ')[0][0] not in ["T","t"]):
                        reading_photons = False; photon_measure = None; photon_unit = None
                    else:
                        active_photon_measure_found = True
                if n == 1:
                    photon_measure_val = float(line.split(' ')[0])
                n += 1
            if reading_electron_source:
                if n == 0:
                    source_fraction = float(line.split(' ')[0])  
                if n == 1:
                    source_energy = float(line.split(' ')[0])    
                if n == 2:
                    source_duration = float(line.split(' ')[0])  
                n+=1            
    print("Time range:",start_t,"-",end_t)
    print("Photon energy:", photon_energy)
    if source_energy is not None:
        print("Source energy:", source_energy)
    param_name_list = ["Energy","Width",photon_measure,"R"]  #TODO Poor format given source energy is now a thing.
    unit_list = [" eV"," fs",photon_unit,""]
    #TODO check that time range is satisfied by files.
    param_dict = dict(
        start_t=start_t,
        end_t=end_t,
        energy=photon_energy,
        width=fwhm,
        fluence=photon_measure_val,
        source_fraction = source_fraction,
        source_energy = source_energy,
        source_duration = source_duration,
    )
    return param_dict, param_name_list,unit_list

def get_sim_elements(handle,input_path=None,molecular_path=None):
    if input_path is None:
        input_path = path.abspath(path.join(__file__ ,"../../input/")) + "/"
    if molecular_path is None:
        molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
        
    molfile = get_mol_file(input_path,molecular_path,handle,"y",out_prefix_text = "Reading elements from")
    outDir = molecular_path + handle 
    intFile = outDir + "/intensity.csv"
    reading = False
    elements = []
    with open(molfile, 'r') as f:
        n = 0
        for line in f:
            if line.startswith("#ATOMS"):
                reading=True            
                continue        
            elif reading == True and (line.startswith("#") or line.startswith("//") or len(line.strip()) == 0):
                break
            if line.startswith("####END####"):
                break
            if reading:
                elem = line.split(' ')[0]      
                elements.append(elem)
    return elements
    

# Contender for world's most convoluted function.
def get_mol_file(input_path, molecular_path, mol, output_mol_query = "",out_prefix_text = "Using"):
    set_highlighted_excepthook()

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
        if y_n_check.casefold() not in map(str.casefold,["n","no"]):  #  If user didn't say to use input...
            if y_n_check.casefold() not in map(str.casefold,["y","yes"]): # If user didn't say to use output either...
                print(unknown_response_qualifier)
            if output_mol_query != "":
                print('\033[95m' + out_prefix_text+" \033[94m'" + os.path.basename(molfile) +"'\033[95m found in output." + '\033[0m')
        else:
            print("Using mol file from " + input_path)
            use_input_mol_file = True
    else: #TODO This case is outdated and useless.
        print("\033[93m[ Missing Mol File ]\033[0m copy of mol file used to generate output not found, searching " +input_path +"\033[0m'" )  
        use_input_mol_file = True
    if use_input_mol_file:       
        molfile = find_mol_file_from_directory(input_path,mol) 
        print(out_prefix_text+" "+ molfile+" found in input.")
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
def get_pdb_paths_dict(my_dir): 
    '''
    my_dir = calling file's directory
    '''
    PDB_PATHS = dict(
        tetra = "targets/5zck.pdb",
        lys = "targets/4et8.pdb", #"targets/2lzm.pdb",
        test = "targets/4et8.pdb", 
        lys_tmp = "targets/4et8.pdb",
        lys_solvated = "solvate_1.0/lys_8_cell.xpdb",
        fcc = "targets/FCC.pdb",
        lys_empty = "targets/lys_points.pdb",
        glycine = "targets/glycine.pdb",
    )      
    for key, value in PDB_PATHS.items():
        PDB_PATHS[key] = my_dir + value
    return PDB_PATHS


def get_spatial_indices(molecular_path, mol):
    output_folder  = molecular_path + mol + '/'
    tmp = "freeDist-"
    indices = []
    for suffix in [f[len(tmp):] for f in os.listdir(output_folder) if f.startswith(tmp)]:
        indices.append(int(suffix.split('.')[0]))
        indices.sort()
    return indices
        

def get_pdb_path(my_dir,key): 
    '''
    my_dir = calling file's directory
    '''
    return get_pdb_paths_dict(my_dir)[key]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_elecs_from_latex(latexlike):
    engine = re.compile(r'(\d[spdf])\^\{(\d+)\}')
    # Parses a chemistry-style configuration to return dict of total number of electrons
    # e.g. $1s^{1}2s^{1}$
    # ignores anything non-numerical or spdf
    qdict = {}
    for match in engine.finditer(latexlike):
        qdict[match.group(1)] = int(match.group(2))
    return qdict

ATOMS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'
       +' Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe').split()
ATOMNO = {}
i = 1
for symbol in ATOMS:
    ATOMNO[symbol] = i
    ATOMNO[symbol + '_fast'] = i
    ATOMNO[symbol + '_faster'] = i
    i += 1
i = 1 
ATOMNO["Gd"] = ATOMNO["Gd_fast"] = ATOMNO["Gd_galli"]  = 64 
ATOMNO["Fe_singleShell"] = ATOMNO["Fe"]
for symbol in list(ATOMNO.keys()):
    if "_" not in symbol:
        ATOMNO[symbol + '_LDA'] = i
        i += 1




def get_data_point(ax,stem,mol_name,mode,SCATTERING_TARGET_DICT,SCATTERING_TARGET,EDGE_SEPARATION,INDEP_VARIABLE,DEP_VARIABLE):
    im_params,_ = SCATTERING_TARGET_DICT[SCATTERING_TARGET]
    indep_variable_key = str(INDEP_VARIABLE)
    dep_variable_key = str(DEP_VARIABLE)
    if INDEP_VARIABLE == EDGE_SEPARATION:
        indep_variable_key += "L" + L_TYPE
    if DEP_VARIABLE is R_FACTOR:
        dep_variable_key = str(DEP_VARIABLE)+"-"+str(SCATTERING_TARGET)
    intensity_averaged = True #TODO
    if DEP_VARIABLE is AVERAGE_CHARGE:
        if intensity_averaged: 
            dep_variable_key += '_avged'
        else:
            dep_variable_key += '_EoP'
    saved_x = saved_y = saved_end_time = None
    dat= get_saved_data(mol_name) 
    if dat is not None:
            saved_x,                              saved_y,                            saved_end_time = (
            dat.get("x").get(indep_variable_key),dat.get("y").get(dep_variable_key),dat.get("end_time")
        )                          
    ## Measure damage 
    if saved_y is None:              
        pl = Plotter(mol_name)
        if mode is AVERAGE_CHARGE:
            intensity_averaged = True
            if intensity_averaged: 
                y = np.average(pl.get_total_charge(atoms="C")*pl.intensityData)/np.average(pl.intensityData)
            else:
                step_index = -1  # Average charge at End-of-pulse 
                y = pl.get_total_charge(atoms="C")[step_index]
        elif mode is R_FACTOR:            
            y = get_R(mol_name,MOLECULAR_PATH,im_params)[0][0]
    else: 
        y = saved_y
    # Get end times (for purposes of checking all simulations ran to completion)
    if saved_end_time is None:
        if saved_y is not None:
            pl = Plotter(mol_name)
        t = pl.timeData[-1] 
    else:
        t = saved_end_time        
    ## Get independent variable (must be last because of sourceless case)
    if INDEP_VARIABLE is PHOTON_ENERGY: 
        energy = get_sim_params(mol_name)[0]["energy"]
    if INDEP_VARIABLE is EDGE_SEPARATION:
        if type(edge_dict[stem][1]) in (list,tuple) and L_TYPE == "one":
            energy = get_sim_params(mol_name)[0]["energy"] - edge_dict[stem][0][2]*1000
        elif type(edge_dict[stem][1]) in (list,tuple) and  L_TYPE == "two":
            energy = get_sim_params(mol_name)[0]["energy"] - edge_dict[stem][0][1]*1000
        elif type(edge_dict[stem][1]) in (list,tuple) and  L_TYPE == "three":
            energy = get_sim_params(mol_name)[0]["energy"] - edge_dict[stem][0][0]*1000
        else:
            energy = get_sim_params(mol_name)[0]["energy"] - max(edge_dict[stem][0])*1000
    if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
        s = get_sim_params(mol_name)[0]
        energy = s["source_energy"]
        if energy is None:
            # No source
            if (s["fluence"],s["width"]) != (fluence,width):
                print(s)
                print((s["fluence"],s["width"]),"!=",(fluence,width))
                raise ValueError("Params do not match expected values for batch")                        
            # Plot horizontal dashed line corresponding to damage with no artificial source.
            ax.axhline(y=y,ls=(5,(10,3)),color=cmap(0),label="No injected")
        elif (s["fluence"],s["width"],s["source_fraction"],s["source_duration"]) != (fluence,width,source_fraction,source_duration):
            print(s)
            print((s["fluence"],s["source_fraction"],s["source_duration"]),"!=",(fluence,width,source_fraction,source_duration))
            raise ValueError("Params do not match expected values for batch")                    
    if energy!= None:
        # Standard case
        x = energy/1000    
    if saved_x is not None:
        assert saved_x == x                                 
    save_dict = {}
    if saved_x is None:
        save_dict["x"] = {indep_variable_key:x}
    if saved_y is None:    
        save_dict["y"] = {dep_variable_key:y}
    if saved_end_time is None:
        save_dict["end_time"] = t
    save_data(mol_name,save_dict,delete_old=True)
    return x,y,t


T_PRECISION = 6 # Truncate past 6 d.p. (millionth of an fs) to Avoid floating point error

class PlotData:
    def __init__(self, abs_molecular_path, mol_name,output_mol_query, max_final_t, max_points,spatial_index,custom_name = None):
        self.molecular_path = abs_molecular_path
        AC4DC_dir = path.abspath(path.join(__file__ ,"../../"))  + "/"
        self.input_path = AC4DC_dir + 'input/'
        molfile = get_mol_file(self.input_path, self.molecular_path, mol_name,output_mol_query,out_prefix_text = "Reading atoms in") 

        # Subplot dictionary
        subplot_name = mol_name.replace('_',' ')
        if custom_name is not None:
            subplot_name = custom_name
        self.target_mol = {'name': subplot_name, 'infile': molfile, 'mtime': path.getmtime(molfile)}

        # Stores the atomic input files read by AC4DC
        self.atomdict = {}
        self.statedict = {}

        # Stores the output of AC4DC
        self.boundData={}
        self.chargeData={}
        self.freeData=None
        self.intensityData=None
        self.energyKnot=None
        self.timeData=None

        self.max_final_t = max_final_t 
        self.max_points = max_points
        self.spatial_tag = "" # This option is for compatibility with simulations from zero-dimensional version of code
        
        self.spatial_index = spatial_index
        if spatial_index is not None:
            self.spatial_tag = f"-{spatial_index}"
        
        self.title_colour = "#4d50b3"
        # Directories of target data
        self.outDir = self.molecular_path + mol_name
        self.freeFile = self.outDir+f"/freeDist{self.spatial_tag}.csv"
        self.intFile = self.outDir + "/intensity.csv"

        self.get_atoms()

        self.raw_int = np.genfromtxt(self.intFile, comments='#', dtype=np.float64) # We reuse this a few times so store it.

        #self.update_outputs()

   # Reads the control file specified by self.mol['infile']
    # and populates the atomdict data structure accordingly
    def get_atoms(self):
        self.atomdict = {}
        with open(self.target_mol['infile'], 'r') as f:
            reading = False
            for line in f:
                if line.startswith("#ATOMS"):
                    reading=True
                    continue
                elif line.startswith("#") or line.startswith("//"):
                    reading=False
                    continue
                if line.startswith("####END####"):
                    break
                if reading:
                    a = line.split(' ')[0].strip()
                    if len(a) != 0:
                        file = self.input_path + 'atoms/' + a + '.inp'
                        self.atomdict[a]={
                            'infile': file,
                            'mtime': path.getmtime(file),
                            'outfile': self.outDir+f"/dist_{a}{self.spatial_tag}.csv"}        
    def get_max_time_range(self):
        return min(self.max_final_t,self.raw_int[-1,0]) - self.raw_int[0,0]
    def set_max_t(self,time_range):
        raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        Q = 10**(-T_PRECISION)
        self.max_final_t = (self.raw_int[0,0] + time_range)//Q*Q  # Truncate to avoid floating point error
        #self.max_final_t = round((self.raw_int[0,0] + time_range)/Q)*Q  # Round to avoid floating point error
    def get_num_usable_points(self):
        return min(len(self.raw_int), self.max_points)

    def update_outputs(self):
        # Get samples of steps separated by the same times.
        times = np.linspace(self.raw_int[0,0],self.max_final_t,self.max_points)
        Q = 10**(-T_PRECISION)
        indices = np.searchsorted(self.raw_int[:,0],(times)//Q*Q)
        #indices = np.searchsorted(self.raw_int[:,0],np.round((times)/Q)*Q )
        indices[1:] -= 1  # index 0 is ignored later anyway.
        np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
        raw = self.raw_int[indices]
        print("Snapshot times set:\n",raw[:,0])
        
        self.intensityData = raw[:,1]
        self.timeData = raw[:, 0]       
        self.energyKnot = np.array(self.get_free_energy_spec(), dtype=np.float64)
        
        raw = np.genfromtxt(self.freeFile, comments='#', dtype=np.float64)
        raw = raw[indices]
        if  len(self.timeData) != len(raw[:,1]):
            raise Exception("time and free lengths don't match")
        self.freeData = raw[:,1:]   
        for a in self.atomdict:
            raw = np.genfromtxt(self.atomdict[a]['outfile'], comments='#', dtype=np.float64)
            self.boundData[a] = raw[:, 1:]
            self.statedict[a] = self.get_bound_config_spec(a)   

    def get_free_energy_spec(self):
        erow = []
        with open(self.freeFile) as f:
            r = csv.reader(f, delimiter=' ')
            for row in r:
                if row[0] != '#':
                    raise Exception('parser could not find energy grid specification, expected #  |')
                reading=False
                for entry in row:
                    # skip whitespace
                    if entry == '' or entry == '#':
                        continue
                    # it's gotta be a pipe or nothin'
                    if reading:
                        erow.append(entry)
                    elif entry == '|':
                        reading=True
                    else:
                        break
                if reading:
                    break
        return erow        
        
    def get_bound_config_spec(self, a):
        # gets the configuration strings corresponding to atom a
        specs = []
        with open(self.atomdict[a]['outfile']) as f:
            r = csv.reader(f, delimiter=' ')
            for row in r:
                if row[0] != '#':
                    raise Exception('parser could not find atomic state specification, expected #  |')
                reading=False
                for entry in row:
                    # skip whitespace
                    if entry == '' or entry == '#':
                        continue
                    # it's gotta be a pipe or nothin'
                    if reading:
                        specs.append(entry)
                    elif entry == '|':
                        reading=True
                    else:
                        break
                if reading:
                    break
        return specs    
    
    def get_density(self, t):
        t_idx = self.timeData.searchsorted(t)
        de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2]) 
        de = de [1:] - de[:-1]
        return np.dot(self.freeData[t_idx, :], de)       

    def aggregate_charges(self,charge_difference=False):
        # populates self.chargeData based on contents of self.boundData
        for a in self.atomdict:
            states = self.statedict[a]
            if len(states) != self.boundData[a].shape[1]:
                msg = 'states parsed from file header disagrees with width of data width'
                msg += ' (got %d, expected %d)' % (len(states), self.boundData[a].shape[1])
                raise RuntimeError(msg)

            initial_charge = 0
            if charge_difference:
                initial_charge = ATOMNO[a] - sum(parse_elecs_from_latex(states[0]).values()) 
            self.chargeData[a] = np.zeros((self.boundData[a].shape[0], ATOMNO[a]+1))
            for i in range(len(states)):
                orboccs = parse_elecs_from_latex(states[i])
                charge = ATOMNO[a] - sum(orboccs.values()) - initial_charge 
                self.chargeData[a][:, charge] += self.boundData[a][:, i]                           
    
    def get_charge(self,atom,charge_difference=False,densities=False,every=1):
        T = self.timeData
        T_start_idx = self.timeData.searchsorted(T[0])
        T_end_idx = self.timeData.searchsorted(T[-1])
        self.aggregate_charges(charge_difference)
        Q = np.zeros(T.shape[0]) # total charge
        colour = None
        assert (atom in self.atomdict)

        # Plot trace for the atom
        atomic_charge = np.zeros(T.shape[0])
        for i in range(self.chargeData[atom].shape[1]):
            atomic_charge += self.chargeData[atom][:,i][T_start_idx:T_end_idx+1]*i

        if not densities:
            atomic_charge /= np.sum(self.chargeData[atom][0]) # average charge

        return atomic_charge

            