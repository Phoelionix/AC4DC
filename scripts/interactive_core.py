import stringprep
import matplotlib.rcsetup as rcsetup
import chart_studio.plotly as py
import numpy as np
# from scipy.interpolate import BSpline
from math import log
import os.path as path
import os
import matplotlib.colors as colors
import sys
# import glob
import csv
import time
import subprocess
import re
from matplotlib.ticker import LogFormatter 
import random
from scipy.optimize import curve_fit
from scipy.stats import linregress
import chart_studio.plotly as py
import plotly.graph_objects as go
import plotly.io as pio
import copy
pio.templates.default = "seaborn" #"plotly_dark" # "plotly"


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


# def get_colors(num, seed):
#     idx = list(np.linspace(0, 1, num))[1:]
#     random.seed(seed)
#     # random.shuffle(idx)
#     idx.insert(0,0)
#     C = plt.get_cmap('nipy_spectral')
#     return C(idx)

class PlotData:
    def __init__(self,target_path, mol_name,output_mol_query):
        self.p = target_path
        molfile = self.get_mol_file(mol_name,output_mol_query) 

        # Subplot dictionary
        subplot_name = mol_name.replace('_',' ')
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

        self.plot_final_t = -9.996#0 # TODO temporary, move to generate_interactive.
        self.max_points = 70 # TODO temporary, move to generate_interactive.
        
        self.title_colour = "#4d50b3"
        # Directories of target data
        self.outDir = self.p + "/output/__Molecular/" + mol_name
        self.freeFile = self.outDir+"/freeDist.csv"
        self.intFile = self.outDir + "/intensity.csv"

        self.get_atoms()

        self.update_outputs()

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
                        
    def update_outputs(self):
        raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        # Remove every other point until we have less than the max number of points in our to-be-plotted time range.
        while self.max_points < np.searchsorted(raw[:,0],self.plot_final_t):
            raw = np.delete(raw, list(range(0, raw.shape[0], 2)),axis=0)       
        self.intensityData = raw[:,1]
        self.timeData = raw[:, 0]       
        self.energyKnot = np.array(self.get_free_energy_spec(), dtype=np.float64)
        
        raw = np.genfromtxt(self.freeFile, comments='#', dtype=np.float64)
        while len(self.timeData) < len(raw[:,1]):
            raw = np.delete(raw, list(range(0, raw.shape[0], 2)),axis=0)
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

    def get_mol_file(self, mol, output_mol_query = ""):
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

class InteractivePlotter:
    # Example initialisation: Plotter(water)
    # --> expects there to be a control file named water.mol within AC4DC/input/ or a subdirectory.
    def __init__(self, target_names, output_mol_query):
        target_path = path.abspath(path.join(__file__ ,"../../")) + "/"
        self.num_plots = len(target_names)
        self.target_data = []
        for mol_name in target_names:    
            self.target_data.append(PlotData(target_path,mol_name,output_mol_query))

        #self.autorun=False  

    # def setup_step_axes(self):
    #     self.fig_steps = py.figure(figsize=(4,3))
    #     self.ax_steps = self.fig_steps.add_subplot(111)
    #     self.fig_steps.subplots_adjust(left=0.18,right=0.97, top=0.97, bottom= 0.16)

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

    def go(self):
        if not self.check_current():
            self.rerun_ac4dc()
            self.update_outputs()
    
    def initialise_interactive(self, target, x_args={}, y_args={}):
        self.fig = go.Figure()
        
        self.fig.update_layout(
            title= target + " - Free-electron distribution",  # Attention: title overwritten by add_time_slider()
            showlegend=False,
            font=dict(
                family="Courier New, monospace",
                size=1,
                # color='rgba(0,0,0,0)' # hide tick numbers? Nope.
            ),
            paper_bgcolor= '#F7CAC9' #'#f7dac9' '#F7CAC9'  '#FFD580' "#92A8D1"  lgrey = '#bbbfbf',
        )
        self.fig.update_xaxes(x_args)
        self.fig.update_yaxes(y_args)        

    # line_args - a list of kwarg dictionaries to be used for the line argument of go.Scatter().
    # colour_mixup - good for distinguishing plots that are on the same timescale.
    def plot_traces(self, saturation = 0.85, normed = True, colour_mixup = True, line_kwargs = [{},{},{},{},{}], fitE = None):
        # Add a group of traces for each target, but the time slider needs to manual separate them out. Probably a better way to do this but oh well.
        for g, target in enumerate(self.target_data):
            min_t = target.timeData[0]   #
            max_t = target.timeData[-1]-(target.timeData[-1]-min_t)/len(target.timeData)  # hack, weird indices fix colours and don't affect data.
            # Add traces, one for each slider step
            X = target.energyKnot
            for j, t in enumerate(target.timeData):
                if t > target.plot_final_t:
                    break

                if j == 0: continue  # Skip empty plot

                data = target.freeData[j,:]
                if normed:
                    tot = target.get_density(t)
                    data /= tot
                    data/=4*3.14

                ### Cute colours  ★ﾟ~(◠ᴗ ◕✿)ﾟ*｡:ﾟ+ 
                rgb = [None,None,None]
                a = 1 
                if len(self.target_data) > 1: a = 0.8
                if colour_mixup:
                    # randomish mix thing
                    # mix = 1- g/len(self.target_data)   
                    # rgb_intensity = [mix*1,0.68 + (1-mix)*0.5,mix*1]  # max = 1
                    # rgb_width = [0.4 + (1-mix)*2, 0.6 + (1-mix)*0.3,0.9 + (1-mix)*0.4]
                    # rgb_bndry = [1,0.6+(1-mix)*0.2,(1-mix)*0.2]

                    # # temp
                    # if g == 1:
                    #     rgb_intensity = [1,0.68,1]
                    #     rgb_width = [0.4 + 2, 0.6,0.9] 
                    #     rgb_bndry = [1,0.6,0]

                    if g == 0:
                    ## blue_grey-mustard
                        rgb_intensity = [0.6,0.6,0.8]  # max = 1
                        rgb_width = [1,2,1]
                        rgb_bndry = [1,1,0]

                        target.title_colour =  "#4d50b3" 
                    
                    else:      
                        target.title_colour =  "#4d50b3"  # "#a44ae8" 
                        ## blue-purp
                        # rgb_intensity = [1,0,1]  # max = 1
                        # rgb_width = [1,0.5,1]
                        # rgb_bndry = [1,0.5,0]

                        ## blue-orange
                        rgb_intensity = [1,0.68,1]  # max = 1
                        rgb_width = [0.5,0.6,0.5]
                        rgb_bndry = [1,0.6,0]                               
                
                else:
                    rgb_intensity = [1,0.68,1]  # max = 1
                    rgb_width = [0.4,0.6,0.9]
                    rgb_bndry = [1,0.6,0]

                # Linear interpolation
                for i in range(len(rgb)):
                    t_norm = (t-min_t)/(max_t-min_t)
                    rgb[i] = saturation * rgb_intensity[i] * (1-((t_norm-rgb_bndry[i])/rgb_width[i])**2)
                    rgb[i] = 255*min(1, max(0, rgb[i]))
                    rgba = tuple(rgb) + (a,)
                col = "rgba" + str(tuple(rgba))      
                
                self.fig.add_trace(
                    go.Scatter(
                        visible=False,
                        # Can't get to work
                        # marker=dict(
                        #     color='LightSkyBlue',
                        #     size=12,
                        # ),
                        # fillpattern = dict(
                        #     shape = "+",                    
                        #     fillmode = 'overlay',
                        # ),
                        line=dict(color=col, **line_kwargs[g]),
                        name="t = " + str(t),
                        x=X,
                        y=data*X))
        
        self.fig.data[0].visible = True

    #----Widgets----#
    # Time Slider
    def add_time_slider(self):
        self.steps_groups = [] # Stores the steps for each plot separately 
        start_step = 0
        time_slider = []
        simul_steps=[]
        for g in range(self.num_plots):
            steps = []
            target = self.target_data[g]
            displaying = "<span style='font-size: 28px; font-family: Roboto'>" +"Displaying:                  </span>"
            subplot_title = dict(text= displaying + "<span style='font-size: 35px;color:"+ target.title_colour +"; font-family: Roboto'>" + target.target_mol["name"]  + "</span>", yanchor = "top", xanchor = "left", pad = dict(b = 0,l=-400))  # margin-top:100px; display:inline-block;
            allplot_title = copy.deepcopy(subplot_title) 
            allplot_title_colour = "#4d50b3" 
            allplot_title["text"] = displaying + "<span style='font-size: 35px;color:"+ allplot_title_colour +"; font-family: Roboto'>" + "          All"
            if g == 0:
                self.fig.update_layout({"title": subplot_title})
            for i in range(len(target.timeData) - 1): # -1 as don't have trace for zeroth time step.
                if target.timeData[i+1] > target.plot_final_t:
                    break                
                step = dict(
                    method="update",
                    args=[
                        {"visible": [False] * len(self.fig.data)},  # style attribute
                        {"title": subplot_title},                   # layout attribute         #Add line below: + '<br>' +  '<span style="font-size: 12px;">line2</span>'}  
                    ],      
                    label= "  " + "%.2f" % target.timeData[i+1],
                )
                step["args"][0]["visible"][start_step + i] = True  #  When at this step toggle i'th trace in target's group to "visible"
                steps.append(step)
                #   Initialise slider that shows all plots.
                if g == 0:
                    simul_step = copy.deepcopy(step)
                    simul_step["label"] = "  " + "%.2f" % target.timeData[i+1]
                    simul_steps.append(simul_step)    
                    simul_step["args"][1]["title"] = allplot_title 
                #   Add later plots' traces at same step (not necessarily same time...).
                elif i < len(simul_steps):
                    simul_steps[i]["args"][0]["visible"][start_step+i] = True    
                    simul_steps[i]["label"] += "  |  " + "%.2f" % target.timeData[i+1]
            ###
            self.steps_groups.append(steps)
            start_step += len(steps)

            time_slider.append(dict(
                active=0,
                tickwidth=0,
                tickcolor = "rgba(0,0,0,0)",
                currentvalue={"prefix": "<span style='font-size: 25px; font-family: Roboto; color = black'>" + target.target_mol["name"] + " - Time [fs]: "},
                pad={"t": 85+90*g,"r": 200,"l":0},
                steps=steps,
                #font = {"color":"rgba(0.5,0.5,0.5,1)"}
            ))
        
        all_slider = dict(
            active = 0,
            tickwidth = 0,
            currentvalue = {"prefix": "<span style='font-size: 25px; font-family: Roboto; color = white;'>" +"All" + " - Times [fs]: "},
            pad={"t": 85+90*(g+1),"r": 200,"l":0},
            steps = simul_steps,
            #font = {"color":"rgba(0.5,0.5,0.5,1)"}
        )

        time_slider.append(all_slider)
        
        self.fig.update_layout(sliders=time_slider)    

    # Scale Button 
    def add_scale_button(self,x_log_args, x_lin_args, y_log_args,y_lin_args,):      
        scale_button = go.layout.Updatemenu(
                buttons=list([
                    dict(
                        args=[{'yaxis': y_log_args, 'xaxis': x_log_args}],
                        label="<br>Log-Log<br>",
                        method="relayout",
                    ),
                    dict(
                        args=[{'yaxis': y_lin_args, 'xaxis': x_log_args}],
                        label="<br>Lin-Log<br>",
                        method="relayout"    
                    ),              
                    dict(
                        args=[{'yaxis': y_lin_args, 'xaxis': x_lin_args}],
                        label="<br>Lin-Lin<br>",
                        method="relayout"              
                    )  
                ]),
                type="buttons",
                direction = "down",
                # anchor top of button to bottom right of graph, then push down.
                pad={"r": 0, "t": 35},
                showactive=True,
                x= 1,
                xanchor="right",
                y= 0, 
                yanchor="top",
                font = {"size": 25,"family": "Roboto"},
        )   
        self.fig.update_layout(
            updatemenus = [scale_button]
        )    

            
    # def plot_fit(self, t, fitE, normed=True, **kwargs):
    #     t_idx = self.timeData.searchsorted(t)
    #     fit = self.energyKnot.searchsorted(fitE)
    #     data = self.freeData[t_idx,:]
    #     if normed:
    #         tot = self.get_density(t)
    #         data /= tot
    #         data/=4*3.14

    #     Xdata = self.energyKnot[:fit]
    #     Ydata = data[:fit]
    #     mask = np.where(Ydata > 0)
    #     T, n = fit_maxwell(Xdata, Ydata)
    #     return self.ax_steps.plot(self.energyKnot, 
    #         maxwell(self.energyKnot, T, n)*self.energyKnot,
    #         '--',label='%3.1f eV' % T, **kwargs)

    # def plot_maxwell(self, kT, n, **kwargs):
    #     return self.ax_steps.plot(self.energyKnot, 
    #         maxwell(self.energyKnot, kT, n)*self.energyKnot,
    #         '--',label='%3.1f eV' % kT, **kwargs)


    # def get_temp(self, t, fitE):
    #     t_idx = self.timeData.searchsorted(t)
    #     fit = self.energyKnot.searchsorted(fitE)
    #     Xdata = self.energyKnot[:fit]
    #     Ydata = self.freeData[t_idx,:fit]
    #     T, n = fit_maxwell(Xdata, Ydata)
    #     return (T, n)

        

# def fit_maxwell(X, Y):
#     guess = [200, 12]
#     # popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess, sigma=1/(X+10))
#     popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess)
#     return popt

# def maxwell(e, kT, n):
#     if kT < 0:
#         return 0 # Dirty silencing of fitting error - note we get negative values from unphysical oscillations, so this increases the average value around this point. -S.P.
#     return n * np.sqrt(e/(np.pi*kT**3)) * np.exp(-e/kT)

# def plot_maxwell(kT, n):
#     e_points = np.logspace(0,4,100)
#     #plt.plot(e_points,maxwell(e_points,kT,n))
#     #pl.ax_steps.plot(e_points,maxwell(e_points,kT,n))
#     pl.ax_steps.plot(e_points,maxwell(e_points,kT,n)*e_points,'--', **kwargs)
#     return

# def lnmaxwell(e, kT, n):
#     return np.log(n) + 0.5*np.log(e/np.pi*kT**3) - e /kT

# def moving_average(a, n=3) :
#     ret = np.cumsum(a, dtype=float)
#     ret[n:] = ret[n:] - ret[:-n]
#     return ret[n - 1:] / n

if __name__ == "__main__":
    pl = Plotter(sys.argv[1])
    # pl.plot_free(log=True,min=1e-7)
    # pl.plot_all_charges()
    plt.show()

