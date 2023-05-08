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
from core_functions import get_mol_file
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
    def __init__(self, abs_molecular_path, mol_name,output_mol_query, max_final_t, max_points,custom_name = None):
        self.molecular_path = abs_molecular_path
        AC4DC_dir = path.abspath(path.join(__file__ ,"../../"))  + "/"
        self.input_path = AC4DC_dir + 'input/'
        molfile = get_mol_file(self.input_path, self.molecular_path, mol_name,output_mol_query) 

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
        
        self.title_colour = "#4d50b3"
        # Directories of target data
        self.outDir = self.molecular_path + mol_name
        self.freeFile = self.outDir+"/freeDist.csv"
        self.intFile = self.outDir + "/intensity.csv"

        self.get_atoms()

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
                            'outfile': self.outDir+"/dist_%s.csv"%a}        
    def get_max_t(self):
        raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        # Get samples of steps separated by the same times. TODO need to fix AC4DC saving points to the nonraw file when loading sim so that the loaded part isnt empty.   
        return min(self.max_final_t,raw[-1,0]) 
    def get_num_usable_points(self):
        raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        return min(len(raw), self.max_points)

    def update_outputs(self):
        raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        # Get samples of steps separated by the same times. TODO need to fix AC4DC saving points to the nonraw file when loading sim so that the loaded part isnt empty.
        times = np.linspace(raw[0,0],self.max_final_t,self.max_points)
        indices = np.searchsorted(raw[:,0],times) 

        np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})
        print("Times plotted:\n",raw[indices,0])
        raw = raw[indices]
        
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

class InteractivePlotter:
    # max_final_t, float, end time in femtoseconds. Not equivalent to time duration
    # max_points, int, number of points (within the timespan) for the interactive to have at maximum.
    def __init__(self, target_names, sim_output_parent_directory, max_final_t = 30, max_points = 70, custom_names = None):
        '''
        output_parent_directory: absolute path
        '''
        self.multi_trace_params = [""]*len(target_names)  # 

        self.num_plots = len(target_names)
        self.target_data = []
        lowest_max_t = np.inf
        lowest_max_points = np.inf
        if custom_names is None:
            custom_names = [None]*self.num_plots
        for i, mol_name in enumerate(target_names):
            custom_name = custom_names[i]
            dat = PlotData(sim_output_parent_directory,mol_name,"y",max_final_t=max_final_t,max_points=max_points,custom_name=custom_name)    
            lowest_max_t = min(dat.get_max_t(),lowest_max_t)
            lowest_max_points = min(lowest_max_points, dat.get_num_usable_points())
            self.target_data.append(dat)
        for dat in self.target_data:
            dat.max_final_t = lowest_max_t
            dat.max_points = lowest_max_points
            dat.update_outputs()

    def set_trace_params(self,target_handle_idx,energy,fwhm,photons):
        self.multi_trace_params[target_handle_idx] = (energy,fwhm,photons)        

    def update_inputs(self):
        self.get_atoms()
        self.mol['mtime'] = path.getmtime(self.mol['infile'])

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
    
    def initialise_interactive(self, plot_title, x_args={}, y_args={}):
        self.fig = go.Figure()
        
        self.fig.update_layout(
            title= plot_title + " - Free-electron distribution",  # Attention: title overwritten by add_time_slider()
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
                if t > target.max_final_t:
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
                    if len(self.target_data) == 2:
                        # special comparison mode...
                        if g == 0:
                        ## blue_grey-mustard
                            rgb_intensity = [0.6,0.6,0.8]  # max = 1
                            rgb_width = [1,2,1]
                            rgb_bndry = [1,1,0]

                            target.title_colour =  "#4d50b3" 
                        
                        else:      
                            target.title_colour =  "#4d50b3"  # "#a44ae8" 
                            ## blue-purp
                            # rgb_intensity = [1,0,1]  
                            # rgb_width = [1,0.5,1]
                            # rgb_bndry = [1,0.5,0]

                            ## blue-orange
                            rgb_intensity = [1,0.68,1]
                            rgb_width = [0.5,0.6,0.5]
                            rgb_bndry = [1,0.6,0]
                    else:
                        #randomish mix thing
                        mix = 1- g/len(self.target_data)
                        rgb_intensity = [mix*1,0.68 + (1-mix)*0.5,mix*1]
                        rgb_width = [0.4 + (1-mix)*2, 0.6 + (1-mix)*0.3,0.9 + (1-mix)*0.4]
                        rgb_bndry = [1,0.6+(1-mix)*0.2,(1-mix)*0.2]
                
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
        simul_step_slider=False
        if self.num_plots > 1:
            # Slider that displays all simulations at the same time step. The limit of plotly without dash.
            simul_step_slider=True
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
                if target.timeData[i+1] > target.max_final_t:
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
                if simul_step_slider:
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
        if simul_step_slider:
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

    # def switching_time_sliders(self):
    #     ## Add the regular time sliders
    #     self.add_time_slider(False)    
    #     # Add buttons to switch sliders
    #     buttons = []
    #     for slider_idx, energy, fwhm, photons  in enumerate(self.multi_trace_params):
    #         vis = False
    #         if slider_idx == 0:
    #             vis = True
    #         button = dict(
    #             label = 'Show first slider',
    #             method = 'update',
    #             args = [
    #                 {'visible': vis},
    #                 {'title': "First Slider", 'xaxis': {'title': 'First X axis'}, 'yaxis': {'title': 'First Y Axis'}, 'sliders': sliders,},
                        
    #             ],
    #             args2 = [
    #                 {'sliders':sliders}
    #             ]
    #         )      
    #         buttons.append(button)    
        #


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


# def lnmaxwell(e, kT, n):
#     return np.log(n) + 0.5*np.log(e/np.pi*kT**3) - e /kT

# def moving_average(a, n=3) :
#     ret = np.cumsum(a, dtype=float)
#     ret[n:] = ret[n:] - ret[:-n]
#     return ret[n - 1:] / n

if __name__ == "__main__":
    raise Exception("No main script")
