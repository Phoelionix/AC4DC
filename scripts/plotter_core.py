import stringprep
import matplotlib.rcsetup as rcsetup
import matplotlib.pyplot as plt
import numpy as np
# from scipy.interpolate import BSpline
from math import log
import os.path as path
import os
import matplotlib.colors as colors
import sys
# import glob
import csv
import subprocess
import re
from matplotlib.ticker import LogFormatter 
import random
from scipy.optimize import curve_fit
from scipy.stats import linregress

#plt.rcParams.update(plt.rcParamsDefault)
#plt.style.use('seaborn-muted')
#plt.style.use('turrell-style')
plt.rcParams["font.size"] = 9

engine = re.compile(r'(\d[spdf])\^\{(\d+)\}')

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_elecs_from_latex(latexlike):
    # Parses a chemistry-style configuration to return total number of electrons
    # e.g. $1s^{1}2s^{1}$
    # ignores anything non-numerical or spdf
    qdict = {}
    for match in engine.finditer(latexlike):
        qdict[match.group(1)] = int(match.group(2))
    return qdict


ATOMS = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'.split()
ATOMNO = {}
i = 1
for symbol in ATOMS:
    ATOMNO[symbol] = i
    ATOMNO[symbol + '_fast'] = i
    i += 1


def get_colors(num, seed):
    idx = list(np.linspace(0, 1, num))[1:]
    random.seed(seed)
    # random.shuffle(idx)
    idx.insert(0,0)
    C = plt.get_cmap('nipy_spectral')
    return C(idx)

class Plotter:
    # Example initialisation: Plotter(water)
    # --> expects there to be a control file named water.mol within AC4DC/input/ or a subdirectory.
    def __init__(self, mol, output_mol_query = ""):
        self.p = path.abspath(path.join(__file__ ,"../../")) + "/"

        molfile = self.get_mol_file(mol,output_mol_query)
        self.mol = {'name': mol, 'infile': molfile, 'mtime': path.getmtime(molfile)}        
        
        # Stores the atomic input files read by ac4dc
        self.atomdict = {}
        self.statedict = {}

        # Outputs
        self.outDir = self.p + "/output/__Molecular/" + mol
        self.freeFile = self.outDir+"/freeDist.csv"
        self.intFile = self.outDir + "/intensity.csv"

        self.boundData={}
        self.chargeData={}
        self.freeData=None
        self.intensityData=None
        self.energyKnot=None
        self.timeData=None

        self.get_atoms()
        self.update_outputs()
        self.autorun=False
        
        self.setup_step_axes()   

    def get_mol_file(self, mol,output_mol_query = ""):
        #### Inputs ####
        #Check if .mol file in outputs
        use_input_mol_file = False
        output_folder  = self.p + 'output/__Molecular/' + mol + '/'
        if not os.path.isdir(output_folder):
            raise Exception("\033[91m Cannot find simulation output folder '" + mol + "'\033[0m" + "(" +output_folder + ")" )
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
                print("Using mol file from " + self.p + 'input/')
                use_input_mol_file = True
        else: 
            print("\033[93m[ Missing Mol File ]\033[0m copy of mol file used to generate output not found, searching input/ directory.\033[0m'" )
            use_input_mol_file = True
        if use_input_mol_file:       
            molfile = self.find_mol_file_from_directory(self.p + 'input/',mol) 
            print("Using: " + molfile)
        return molfile
    
    def find_mol_file_from_directory(self, input_directory, mol):
        # Get molfile from all subdirectories in input folder.  Old comment -> #molfile = self.p+"/input/"+mol+".mol"
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

    def setup_step_axes(self):
        self.fig_steps = plt.figure(figsize=(4,3))
        self.ax_steps = self.fig_steps.add_subplot(111)
        self.fig_steps.subplots_adjust(left=0.18,right=0.97, top=0.97, bottom= 0.16)
        

    # Reads the control file specified by self.mol['infile']
    # and populates the atomdict data structure accordingly
    def get_atoms(self):
        self.atomdict = {}
        with open(self.mol['infile'], 'r') as f:
            reading = False
            for line in f:
                if line.startswith("#ATOMS"):
                    reading=True
                    continue
                elif line.startswith("#"):
                    reading=False
                    continue
                if reading:
                    a = line.split(' ')[0].strip()
                    if len(a) != 0:
                        file = self.p + '/input/atoms/' + a + '.inp'
                        self.atomdict[a]={
                            'infile': file,
                            'mtime': path.getmtime(file),
                            'outfile': self.outDir+"/dist_%s.csv"%a}

    def update_inputs(self):
        self.get_atoms()
        self.mol['mtime'] = path.getmtime(self.mol['infile'])

    def rerun_ac4dc(self):
        cmd = self.p+'/bin/ac4dc2 '+self.mol['infile']
        print("Running: ", cmd)
        subprocess.run(cmd, shell=True, check=True)

    def check_current(self):
        # Pull most recent atom mod time
        self.update_inputs()
        # Outputs are made at roughly the same time: Take the oldest.
        out_time = path.getmtime(self.freeFile)
        for atomic in self.atomdict.values():
            out_time = min(out_time, path.getmtime(atomic['outfile']))

        if self.mol['mtime'] > out_time:
            print("Master file %s is newer than most recent run" % self.mol['infile'])
            return False

        for atomic in self.atomdict.values():
            if atomic['mtime'] > out_time:
                print("Dependent input file %s is newer than most recent run" % atomic['infile'])
                return False
        return True

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


    def get_atomic_numbers(self):
        # Get atomic number by considering states at time 0.
        return_dict = {}
        for a in self.atomdict:
            states = self.statedict[a]  
            atomic_number = 0
            for val in parse_elecs_from_latex(states[0]).values():
                atomic_number += val
            return_dict[a] = atomic_number
        return return_dict
    


    def initialise_coherence_params(self, start_t,end_t, q_max, photon_energy, q_fineness=50,t_fineness=50,naive=False):
        args = locals()
        self.__dict__ .update(args)  # Is this a coding sin?
        self.q_max *= 2*np.pi  # Change to wavenumber.  TODO change code so that the form factor calcs are the only part that use this convention.
 
    
    # q in units of bohr^-1
    def theta_to_q(self,scattering_angle_deg,photon_energy=None):
        if photon_energy == None:
            photon_energy = self.photon_energy
        theta = scattering_angle_deg*np.pi/180
        bohr_per_nm = 18.8973; c_nm = 2.99792458e17; h = 6.582119569e-16 *2 *np.pi
        wavelength =  h*c_nm/photon_energy
        wavelength *= bohr_per_nm
        return 4*np.pi*np.sin(theta)/wavelength
    # q in units of bohr^-1
    def q_to_theta(self,q,photon_energy=None):
        if photon_energy == None:
            photon_energy = self.photon_energy        
        bohr_per_nm = 18.8973; c_nm = 2.99792458e17; h = 6.582119569e-16 *2 *np.pi
        wavelength =  h*c_nm/photon_energy   
        wavelength *= bohr_per_nm
        return 180/np.pi*np.arcsin(q*wavelength/(4*np.pi))     


    def get_A_bar_matrix(self,atom1,atom2):
        # Matrix constructor to be vectorised and passed to np.fromfunction
        def A_maker(i,j,atom1,atom2):
            # Transform indices to desired q
            q1 = i*self.q_max/self.q_fineness
            q2 = j*self.q_max/self.q_fineness
            if self.naive == False:
                return self.get_A_bar(self.start_t,self.end_t,q1,q2,atom1,atom2,self.t_fineness)
            elif self.naive == True:
                return self.get_naive_A_bar(self.start_t,self.end_t,q1,q2,atom1,atom2,self.t_fineness)         
                   
        A_norm = 1/A_maker(0,0,atom1,atom2) 
        return A_norm*np.fromfunction(np.vectorize(A_maker),(self.q_fineness,self.q_fineness,),atom1=atom1,atom2=atom2)
    
    def plot_R_factor(self, atoms_subset=None):
        atoms = atoms_subset
        if atoms_subset == None:
            atoms = self.atomdict

        I_real = None
        I_ideal = None

        ### I_real
        for atom1 in atoms:
            for atom2 in atoms:
                A_matrix = self.get_A_bar_matrix(atom1,atom2)
                if I_real == None:
                    I_real = A_matrix
                else:
                    I_real += A_matrix
        ### I_ideal.
        # Using naive setting (so that np.trapz doesn't return 0), and start_time, A_bar will just be the ideal factor squared
        # Temporarily overwrite variables
        old_naive_param = self.naive
        old_end_t_param = self.end_t
        self.naive = True
        self.end_t = self.start_t
        for atom1 in atoms:
            for atom2 in atoms:
                A_matrix_ideal = self.get_A_bar_matrix(atom1,atom2)
                if I_ideal == None:
                    I_ideal = A_matrix_ideal
                else:
                    I_ideal += A_matrix_ideal
        self.naive = old_naive_param
        self.end_t = old_end_t_param
        ### Weird R_factor_inspired_matrix thing that probably isn't a thing because I got confused by a paper's notation:
        R_like_matrix_that_is_not_a_thing =  (np.sqrt(I_real) - np.sqrt(I_ideal))/np.sqrt(I_ideal)
        
        # Plot R over different q ranges, up to q_max which corresponds to the best resolution but should have the lowest R-factor.
        
        
        q_lengths = np.linspace(0,self.q_max,self.q_fineness-1)  #TODO make better and don't try and do in head...     

        R_factor = []
        # Compute single term of R's sum
        def get_R_num_term(idx):
            I_r = I_real[idx][idx]
            I_i = I_ideal[idx][idx]
            return abs(np.sqrt(I_r)-np.sqrt(I_i))
        def get_R_den_term(idx):
            I_i = I_ideal[idx][idx]    
            return np.sqrt(I_i)        
        # Get R for each q_length
        cumulative_R_num = 0
        cumulative_R_den = 0
        for i,q in enumerate(q_lengths): 
            cumulative_R_num += get_R_num_term(i)
            cumulative_R_den += get_R_den_term(i)
            R_factor.append(cumulative_R_num/cumulative_R_den)  
        
        q_resolution = np.array([[3.2,1],[0.032,100]])  # Based off Hau-Riege 2007 TODO research relationship between q and resolution properly.
        k_lengths = q_lengths/2/np.pi    #TODO This is a mess. really need to deal with this k and q thing properly.
        q_theta = []
        for i, k in enumerate(k_lengths):
            if i%(int(self.q_fineness/5)) == 1:
                q_theta.append([k,self.q_to_theta(k)])
        q_theta.reverse()
        q_theta = np.array(q_theta)
        
        dual_axis = q_theta
        #top_label = r"Resolution (\AA)"
        top_label = r"Theta"         
        
        fig_size = 5
        self.fig = plt.figure(figsize=(fig_size,fig_size))
        ax = self.fig.add_subplot(111)
        ax.set_xlabel("k length",fontsize = 12)
        ax.set_xscale("log")
        ax.set_xlim(dual_axis[0][0],dual_axis[-1][0])
        ax.set_ylim(0,0.1) 
        ax.set_ylabel("R")
        ax.plot(k_lengths,R_factor)

        # Resolution or angle axis 
        top_tick_locations = dual_axis[:,0]
        top_tick_labels = dual_axis[:,1]    
        top_tick_labels = ['%.2f' % l for l in top_tick_labels]
         
        ax2 = ax.twiny()              
        #ax2.set_xlim(q_resolution[1][1],q_resolution[0][1])
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticklabels(top_tick_labels)
        ax2.set_xticks(top_tick_locations)
        ax2.set_xlabel(top_label)
        
    

    def plot_A_map(self, atom1, atom2, vmin=0.,vmax=1.,title=r"Foreground coherence $\bar{A}$"):
        q_fineness = self.q_fineness
        
        A_matrix = self.get_A_bar_matrix(atom1,atom2)

        fig_size = 5
        self.fig = plt.figure(figsize=(fig_size,fig_size))
        ax = self.fig.add_subplot(111)
        ax.set_title(title)
        ax.set_xlabel("$k_{1}$",fontsize=12)
        ax.set_ylabel("$k_{2}$", rotation=0,labelpad = 7,fontsize=12)
        ax.set_xlim(0,q_fineness-1)  
        ax.set_ylim(0,q_fineness-1)  
        ticks = np.linspace(0,q_fineness-1,5)
        ticklabels = ["{:6.2f}".format(i) for i in ticks/(q_fineness-1)*self.q_max/(2*np.pi)]
        
        
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticklabels)  

        im = ax.imshow(
            A_matrix,
            cmap = plt.get_cmap('PiYG'),
            vmin=vmin,
            vmax=vmax,
        )
        cbar = self.fig.colorbar(im, ax=ax, extend='both')
        cbar.minorticks_on()
        

    # Get the foreground coherence matrix A bar.(Martin A.V. & Quiney H.M., 2016)
    # a1, str, atomic name key
    # a2, str, other atomic name key
    def get_A_bar(self, start_t, end_t, q1, q2, atom1, atom2, t_fineness = 50):
        # f_1 = scattering factor of atomic species 1.
        #Integrand = <(f_1)(f_2)*> = (f_Z1)(f_Z2) since real scatt. factors   TODO check where imaginary f may be relevant
        
        def integrand(idx):
            # We pass in indices from 0 to fineness-1, transform to time:
            t = start_t + idx/t_fineness*(end_t-start_t) 
            val = self.get_form_factor(q1,[atom1],t)*self.get_form_factor(q2,[atom2],t)
            return val
        form_factor_product = np.fromfunction(integrand,(t_fineness,))   # Happily it works without using np.vectorise, which is far more costly.
        time = np.linspace(start_t,end_t,t_fineness)
        #Approximate integral with composite trapezoidal rule.
        A_bar = np.trapz(form_factor_product,time)
        return A_bar

    def get_naive_A_bar(self, start_t, end_t, q1, q2, atom1, atom2, t_fineness = 50):
        time = np.linspace(start_t,end_t,t_fineness)
        form_factor_1 = self.get_form_factor(q1,[atom1],time)  
        form_factor_2 = self.get_form_factor(q2,[atom2],time)

        #Average form factors 
        average_form_factor_1 = np.average(form_factor_1)
        average_form_factor_2 = np.average(form_factor_2)
        A_bar = average_form_factor_1 * average_form_factor_2
        return A_bar        

    def plot_form_factor(self,times = [-7.5],q_max = None, atoms=None):  # q = 4*np.pi = 2*k. c.f. Sanders plot.
        if q_max == None:  #TODO check if can remove argument entirely.
            q_max = self.q_max
        if atoms == None:
            atoms = self.atomdict  # All atoms

        #lattice_const = 18.11**(1/3)


        #TODO figure out overlaying.
        q = np.linspace(0,q_max,250)
        k = q/2/np.pi  # Convert to momentum TODO assuming this is also units of lattice parameter, but need to check this is right, seems to match...
        k_max = q_max /2/np.pi
        self.fig = plt.figure(figsize=(3,2.5))
        ax = self.fig.add_subplot(111)

        cmap=plt.get_cmap('plasma')
        min_time = np.inf
        max_time = -np.inf 
        for time in times:
            min_time = min(min_time,time)
            max_time = max(max_time,time) 

        f_avg =[]
        for time in times:
            f = [self.get_form_factor(x,atoms,time=time) for x in q]
            f_avg.append(f)
            ax.plot(k, f,label="t = " + str(time)+" fs",color=cmap((time-min_time)/(0.0001+max_time-min_time)))     
        f_avg = np.average(f_avg,axis=0)
        ax.plot(k, f_avg,'k--',label="Average")  

        ax.set_title("")
        ax.set_xlabel("k")
        ax.set_ylabel("Form factor (normed)")
        ax.set_xlim(0,k_max)   
        ax.legend()
        self.fig.set_size_inches(6,5)           

        #core_f = [self.get_form_factor(x,atomic_numbers,time=time,n=1) for x in q]
        #ax.plot(q,core_f)
    
    def get_form_factor(self,q,atoms,time=-7.5,n=None):
        idx = np.searchsorted(self.timeData,time)
         
        # Iterate through each a passed. 
        for a in atoms:
            ff = 0
            states = self.statedict[a]   
            atomic_density = self.boundData[a][0, 0]
            for i in range(len(states)):
                if n != None:
                    if n != i+1:
                        continue
                orboccs = parse_elecs_from_latex(states[i])             
                state_density = self.boundData[a][idx, i]       
                # Get the form factors for each shell
                occ_list = [-99]*10
                for orb, occ in orboccs.items():
                    l = int(orb[0]) - 1
                    if occ_list[l] == -99:
                        occ_list[l] = 0
                    occ_list[l] += occ
                occ_list = occ_list[:len(occ_list)-occ_list.count(-99)]
                shielding = SlaterShielding(self.atomic_numbers[a],occ_list)              
                ff += shielding.get_atomic_ff(q,state_density,atomic_density)
        return ff
                # Get the average atomic form factor


    def print_bound_slice(self,time=-7.5):
        idx = np.searchsorted(self.timeData,time)
        for a in self.atomdict:
            states = self.statedict[a]

            print("Atom:",a)

            atomic_number = 0

            # Get atomic number and density by considering states at time 0.
            
            for val in parse_elecs_from_latex(states[0]).values():
                atomic_number += val
            atomic_density = self.boundData[a][0, 0]
            print("atomic number, density",atomic_number,atomic_density)
            
            average_occupancy = dict.fromkeys(parse_elecs_from_latex(states[0]), 0)
            for i in range(len(states)):
                orboccs = parse_elecs_from_latex(states[i])              
                state_density = self.boundData[a][idx, i]
                print("Occs: ", orboccs ,", time: ",self.timeData[idx], ", density: ",state_density)
                # Get the average number of orbitals in each state
                for orb, occ in orboccs.items(): 
                    average_occupancy[orb] += occ*(state_density/atomic_density)

            print("Average occupancy:",average_occupancy)
            print("Atomic density:", atomic_density)

    def aggregate_charges(self):
        # populates self.chargeData based on contents of self.boundData
        for a in self.atomdict:
            states = self.statedict[a]
            if len(states) != self.boundData[a].shape[1]:
                msg = 'states parsed from file header disagrees with width of data width'
                msg += ' (got %d, expected %d)' % (len(states), self.boundData[a].shape[1])
                raise RuntimeError(msg)

            self.chargeData[a] = np.zeros((self.boundData[a].shape[0], ATOMNO[a]+1))
            for i in range(len(states)):
                orboccs = parse_elecs_from_latex(states[i])
                charge = ATOMNO[a] - sum(orboccs.values())
                self.chargeData[a][:, charge] += self.boundData[a][:, i]


    def update_outputs(self):
        raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        self.intensityData = raw[:,1]
        self.timeData = raw[:, 0]
        self.energyKnot = np.array(self.get_free_energy_spec(), dtype=np.float64)
        raw = np.genfromtxt(self.freeFile, comments='#', dtype=np.float64)
        self.freeData = raw[:,1:]
        for a in self.atomdict:
            raw = np.genfromtxt(self.atomdict[a]['outfile'], comments='#', dtype=np.float64)
            self.boundData[a] = raw[:, 1:]
            self.statedict[a] = self.get_bound_config_spec(a)
        self.atomic_numbers = self.get_atomic_numbers()

    def go(self):
        if not self.check_current():
            self.rerun_ac4dc()
            self.update_outputs()

    # makes a blank plot showing the intensity curve
    def setup_axes(self):
        self.fig = plt.figure(figsize=(3,2.5))
        ax = self.fig.add_subplot(111)
        ax2 = ax.twinx()
        ax2.plot(self.timeData, self.intensityData, lw = 1, c = 'black', ls = ':', alpha = 0.7)
        # ax2.set_ylabel('Pulse Intensity (photons cm$^{-2}$ s$^{-1}$)')
        ax.set_xlabel("Time (fs)")
        ax.tick_params(direction='in')
        ax2.axes.get_yaxis().set_visible(False)
        # ax2.tick_params(None)
        ax.get_xaxis().get_major_formatter().labelOnlyBase = False
        return (ax, ax2)

    def plot_atom_total(self, a):
        ax, _ax2 = self.setup_axes()
        tot = np.sum(self.boundData[a], axis=1)
        ax.plot(self.timeData, tot)
        ax.set_title("Configurational dynamics")
        ax.set_ylabel("Density")
        self.fig.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)
        
    def plot_atom_raw(self, a):
        ax, _ax2 = self.setup_axes()
        for i in range(self.boundData[a].shape[1]):
            ax.plot(self.timeData, self.boundData[a][:,i], label = self.statedict[a][i])
        ax.set_title("Configurational dynamics")
        ax.set_ylabel("Density")
        self.fig.subplots_adjust(left=0.2, right=0.92, top=0.93, bottom=0.1)

    # Seems to plot ffactors through time, "timedata" is where the densities come from - the form factor text file corresponds to specific configurations (vertical axis) at different q (horizontal axis)
    # In any case, my ffactor function is probably bugged since it doesn't match this. - S.P.
    def plot_ffactor_time_slices(self, a, num_tsteps = 10, timespan = None, show_avg = True, **kwargs):

        if timespan is None:
            timespan = (self.timeData[0], self.timeData[-1])

        start_idx = self.timeData.searchsorted(timespan[0])
        stop_idx = self.timeData.searchsorted(timespan[1])

        ff_path = path.abspath(path.join(__file__ ,"../../output/"+a+"/Xsections/Form_Factor.txt"))
        fdists = np.genfromtxt(ff_path)
        # These correspond to the meaning of the FormFactor.txt entries themselves
        KMIN = 0
        KMAX = 2
        dim = len(fdists.shape)
        kgrid = np.linspace(KMIN,KMAX,fdists.shape[0 if dim == 1 else 1])
        fig2 = plt.figure()
        ax = fig2.add_subplot(111)
        ax.set_xlabel('$k$ (atomic units)')
        ax.set_ylabel('Form factor (arb. units)')
        
        timedata = self.boundData[a][:,:-1] # -1 excludes the bare nucleus
        temporary_fact = 3/0.994/0.99992  # Just normalising carbon - S.P.
        dynamic_k = temporary_fact*np.tensordot(fdists.T, timedata.T,axes=1) 
        step = (stop_idx - start_idx) // num_tsteps
        cmap=plt.get_cmap('plasma')
        fbar = np.zeros_like(dynamic_k[:,0])

        n=0
        for i in range(start_idx, stop_idx, step):
            ax.plot(kgrid, dynamic_k[:,i], label='%1.1f fs' % self.timeData[i], color=cmap((i-start_idx)/(stop_idx - start_idx)))
            fbar += dynamic_k[:,i]
            n += 1

        fbar /= n
        if show_avg:
            ax.plot(kgrid, fbar, 'k--', label=r'Effective Form Factor')
        freal = dynamic_k[:,0]
        print("R = ", np.sum(np.abs(fbar - freal))/np.sum(freal))
        freal /= np.sum(freal)
        fbar /= np.sum(fbar)
        print("Normed R = ", np.sum(np.abs(fbar - freal))/np.sum(freal))
        return (fig2, ax)

    def plot_charges(self, ax, a, rseed=404):
        self.aggregate_charges()

        #print_idx = np.searchsorted(self.timeData,-7.5)

        ax.set_prop_cycle(rcsetup.cycler('color', get_colors(self.chargeData[a].shape[1],rseed)))
        for i in range(self.chargeData[a].shape[1]):
            max_at_zero = np.max(self.chargeData[a][0,:])
            mask = self.chargeData[a][:,i] > max_at_zero*2
            mask |= self.chargeData[a][:,i] < -max_at_zero*2
            Y = np.ma.masked_where(mask, self.chargeData[a][:,i])
            ax.plot(self.timeData, Y, label = "%d+" % i)
            #print("Charge: ", i ,", time: ",self.timeData[print_idx], ", density: ",self.chargeData[a][print_idx,i])
        # ax.set_title("Charge state dynamics")
        ax.set_ylabel(r"Density (\AA$^{-3}$)")

        self.fig.legend(loc = "right")
        self.fig.subplots_adjust(left=0.11, right=0.81, top=0.93, bottom=0.1)

    def plot_subshell(self, a, subshell='1s',rseed=404):
        if not hasattr(self, 'ax_subshell'):
            self.ax_subshell, _ax2 = self.setup_axes()

        states = self.statedict[a]
        tot = np.zeros_like(self.boundData[a][:,0])
        for i in range(len(states)):
            orboccs = parse_elecs_from_latex(states[i])
            tot += orboccs[subshell]*self.boundData[a][:,i]

        self.ax_subshell.plot(self.timeData, tot, label=subshell)
        # ax.set_title("Charge state dynamics")
        self.ax_subshell.set_ylabel(r"Electron density (\AA$^{-3}$)")

        self.fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2)


    def plot_tot_charge(self, every=1):
        ax, _ax2 = self.setup_axes()
        self.aggregate_charges()
        self.fig.subplots_adjust(left=0.22, right=0.95, top=0.95, bottom=0.17)

        T = self.timeData[::every]
        self.Q = np.zeros(T.shape[0])
        for a in self.atomdict:
            atomic_charge = np.zeros(T.shape[0])
            for i in range(self.chargeData[a].shape[1]):
                atomic_charge += self.chargeData[a][::every,i]*i
            ax.plot(T, atomic_charge, label = a)
            self.Q += atomic_charge

        de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])
        de = de [1:] - de[:-1]
        tot_free_Q =-1*np.dot(self.freeData, de)
        ax.plot(T, tot_free_Q[::every], label = 'Free')
        ax.set_ylabel("Charge density ($e$ \AA$^{-3}$)")
        self.Q += tot_free_Q[::every]
        # ax.set_title("Charge Conservation")
        ax.plot(T, self.Q, label='total')
        ax.legend(loc = 'upper left')
        return ax


    def plot_all_charges(self, rseed=404):
        ax, _ax2 = self.setup_axes()
        for a in self.atomdict:
            self.plot_charges(ax, a, rseed)

    def plot_free(self, N=100, log=False, min = 0, max=None, every = None):
        self.fig_free = plt.figure(figsize=(3.1,2.5))
        ax = self.fig_free.add_subplot(111)
        self.fig_free.subplots_adjust(left=0.12, top=0.96, bottom=0.16,right=0.95)

        if every is None:
            Z = self.freeData.T
            T = self.timeData
        else:
            Z = self.freeData.T [:, ::every]
            T = self.timeData [::every]
        
        Z = np.ma.masked_where(Z < min, Z)

        if max is not None:
            Z = np.ma.masked_where(Z > max, Z)
        else:
            max = Z.max()

        # if log:
        #     Z = np.log10(Z)

        
        ax.set_facecolor('black')

        if log:
            cm = ax.pcolormesh(T, self.energyKnot*1e-3, Z, shading='gouraud',norm=colors.LogNorm(vmin=min, vmax=max),cmap='magma',rasterized=True)
            cbar = self.fig_free.colorbar(cm)
        else:
            cm = ax.contourf(T, self.energyKnot, Z, N, cmap='magma',rasterized=True)
            cbar = self.fig_free.colorbar(cm)

        ax.set_ylabel("Energy (keV)")
        ax.set_xlabel("Time (fs)")
        
        # if log:
        #     minval = np.floor(np.min(Z, axis=(0,1)))
        #     maxval = np.ceil(np.max(Z,axis=(0,1)))
        #     formatter = LogFormatter(10, labelOnlyBase=True) 
            
        #     cbar.ax.yaxis.set_ticks(vals)
        #     cbar.ax.yaxis.set_ticklabels(10**vals)
        # else:
        #     cbar = self.fig_free.colorbar(cm)
        cbar.ax.set_ylabel('Free Electron Density, Å$^{-3}$', rotation=270,labelpad=20)

        # plot the intensity
        ax2 = ax.twinx()
        ax2.plot(T, self.intensityData[::every],'w:',linewidth=1)

        ax2.get_yaxis().set_visible(False)

        

    def plot_free_raw(self, N=100, log=False, min = None, max=None):
        plt.figure()
        rawdata = np.genfromtxt(self.outDir+'/freeDistRaw.csv')
        T = rawdata[:,0]
        self.rawZ = rawdata[:,1:].T
        Z = self.rawZ

        if log:
            Z = np.log(Z)
        
        if min is not None:
            Z = np.ma.masked_where(Z < min, Z)

        if max is not None:
            Z = np.ma.masked_where(Z > max, Z)

        Yax = np.arange(Z.shape[0])
        plt.contourf(T, Yax, Z, N, shading='nearest',cmap='magma')
        plt.title("Free electron energy distribution")
        plt.ylabel("Energy (eV)")
        plt.xlabel("Time, fs")
        plt.show()
        plt.colorbar()

    # Plots a single point in time.
    def plot_step(self, t, normed=True, fitE=None, **kwargs):        
        self.ax_steps.set_xlabel('Energy (eV)')
        self.ax_steps.set_ylabel('$f(\\epsilon) \\Delta \\epsilon$')
        self.ax_steps.loglog()
        # norm = np.sum(self.freeData[n,:])
        n = self.timeData.searchsorted(t)
        data = self.freeData[n,:]
        X = self.energyKnot

        if normed:
            tot = self.get_density(t)
            data /= tot
            data/=4*3.14
        
        return self.ax_steps.plot(X, data*X, label='%1.1f fs' % t, **kwargs)

    def plot_fit(self, t, fitE, normed=True, **kwargs):
        t_idx = self.timeData.searchsorted(t)
        fit = self.energyKnot.searchsorted(fitE)
        data = self.freeData[t_idx,:]
        if normed:
            tot = self.get_density(t)
            data /= tot
            data/=4*3.14

        Xdata = self.energyKnot[:fit]
        Ydata = data[:fit]
        mask = np.where(Ydata > 0)
        T, n = fit_maxwell(Xdata, Ydata)
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, T, n)*self.energyKnot,
            '--',label='%3.1f eV' % T, **kwargs)

    def plot_maxwell(self, kT, n, **kwargs):
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, kT, n)*self.energyKnot,
            '--',label='%3.1f eV' % kT, **kwargs)


    def get_temp(self, t, fitE):
        t_idx = self.timeData.searchsorted(t)
        fit = self.energyKnot.searchsorted(fitE)
        Xdata = self.energyKnot[:fit]
        Ydata = self.freeData[t_idx,:fit]
        T, n = fit_maxwell(Xdata, Ydata)
        return (T, n)

    def get_density(self, t):
        t_idx = self.timeData.searchsorted(t)
        de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])
        de = de [1:] - de[:-1]
        return np.dot(self.freeData[t_idx, :], de)

        

def fit_maxwell(X, Y):
    guess = [200, 12]
    # popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess, sigma=1/(X+10))
    popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess)
    return popt

def maxwell(e, kT, n):
    if kT < 0:
        return 0 # Dirty silencing of fitting error - note we get negative values from unphysical oscillations, so this increases the average value around this point. -S.P.
    return n * np.sqrt(e/(np.pi*kT**3)) * np.exp(-e/kT)

def plot_maxwell(kT, n):
    e_points = np.logspace(0,4,100)
    #plt.plot(e_points,maxwell(e_points,kT,n))
    #pl.ax_steps.plot(e_points,maxwell(e_points,kT,n))
    pl.ax_steps.plot(e_points,maxwell(e_points,kT,n)*e_points,'--', **kwargs)
    return

def lnmaxwell(e, kT, n):
    return np.log(n) + 0.5*np.log(e/np.pi*kT**3) - e /kT

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Slater form factor
# Warning: Only for states with exclusively s, p orbitals
class SlaterShielding:
    # Z:  nuclear charge
    # shell_config: The number of electrons in each shell, passed as a list e.g. [1,8,8]
    def __init__(self,Z,shell_occs):
        self.Z = Z 
        self.shell_occs = shell_occs
        self.s_p_slater_shielding()
    def s_p_slater_shielding(self):
        self.s = []
        for n in range(len(self.shell_occs)):  
            s = 0
            for m, N_e in enumerate(self.shell_occs): # N_e = num (screening) electrons
                if m > n:
                    continue
                if m == n:
                    s_factor = 0.35
                    if n == 0:
                        s_factor = 0.30
                    s += s_factor * max(0, N_e - 1)
                elif m == n - 1:
                    s += 0.85*N_e
                else:
                    s += N_e           
            self.s.append(s) 
            if s > self.Z:
                print("Warning, s =",self.s)        
                    
    def get_shell_ff(self, q, shell_num ):
        lamb = self.Z - self.s[shell_num-1]
        lamb/=shell_num
        D = (4*lamb**2+q**2)
        if shell_num == 1:            
            return 16*lamb**4/D**2
        if shell_num == 2:
            return 64*lamb**6*(4*lamb**2-q**2)/D**4
        if shell_num == 3:
            return 128*lamb**7/6*(192*lamb**5-160*lamb**3*q**2+12*lamb*q**4)/D**6
        else:
            print("ERROR! shell_ff")
        
    # Gets a single configuration's form factor * its density.
    def get_atomic_ff(self, q, density, atomic_density):
        ff = 0
        norm = 1/atomic_density/self.Z # such that at t = 0, form factor = 1. (At t = 0, we have neutral config density = atomic density, and a factor of N where N is the number of shells, since the shell ff is normalised.)
        for i in range(len(self.shell_occs)):
            ff += self.get_shell_ff(q,i+1)*norm*self.shell_occs[i]
        return ff * density

if __name__ == "__main__":
    pl = Plotter(sys.argv[1])
    # pl.plot_free(log=True,min=1e-7)
    # pl.plot_all_charges()
    plt.show()
