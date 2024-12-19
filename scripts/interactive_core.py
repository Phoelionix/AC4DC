import stringprep
import matplotlib.rcsetup as rcsetup
import chart_studio.plotly as py
import numpy as np
# from scipy.interpolate import BSpline
from math import log
import os.path as path
import matplotlib.colors as colors
# import glob
import csv
from matplotlib.ticker import LogFormatter 
from scipy.optimize import curve_fit
from scipy.stats import linregress
import chart_studio.plotly as py
import plotly.graph_objects as go
import plotly.io as pio
import copy
from core_functions import PlotData, get_mol_file, parse_elecs_from_latex, get_spatial_indices,ATOMS, ATOMNO, T_PRECISION
import matplotlib as plt
pio.templates.default = "seaborn" #"plotly_dark" # "plotly"


# ATOMS = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'.split()
# ATOMNO = {}
# i = 1
# for symbol in ATOMS:
#     ATOMNO[symbol] = i
#     ATOMNO[symbol + '_fast'] = i
#     #ATOMNO[symbol + '_faster'] = i
#     i += 1
# ATOMNO["Gd"] = i 
# ATOMNO["Gd_galli"] = i 
# ATOMNO["Gd_fast"] = i
# i+= 1


# def get_colors(num, seed):
#     idx = list(np.linspace(0, 1, num))[1:]
#     random.seed(seed)
#     # random.shuffle(idx)
#     idx.insert(0,0)
#     C = plt.get_cmap('nipy_spectral')
#     return C(idx)


class InteractivePlotter:
    # max_final_t, float, end time in femtoseconds. Not equivalent to time duration
    # max_points, int, number of points (within the timespan) for the interactive to have at maximum.
    def __init__(self, target_names, sim_output_parent_directory, max_final_t = 30, max_points = 70, spatial_indices=None, custom_names = None,hide_legend=False,legend_title=None,use_electron_density = False,presentation_mode=False,font_size=35,legend_font_size=35,times_in_legend=True,inset = False):
        '''
        output_parent_directory: absolute path
        use_electron_density: If True, plot electron density rather than energy density
        Various changes, such as bigger font, etc. for presentation purposes.
        '''
        self.multi_trace_params = [""]*len(target_names)  # 
        self.use_electron_density = use_electron_density
        self.presentation_mode = False
        if presentation_mode:
            self.presentation_mode = True

        self.font_size = font_size
        self.legend_font_size = legend_font_size
        self.hide_legend = hide_legend
        self.times_in_legend = times_in_legend
        self.legend_title = legend_title
        self.inset = inset
        if spatial_indices is None:
            spatial_indices = get_spatial_indices( sim_output_parent_directory, target_names[0])
            for target in target_names:
                assert(spatial_indices==get_spatial_indices( sim_output_parent_directory, target))         
        
        self.num_plots = len(target_names)*len(spatial_indices)
        if custom_names is None:
            custom_names = [None]*self.num_plots
        self.input_data_args = {
            "target_names":target_names,
            "custom_names":custom_names,
            "sim_output_parent_directory":sim_output_parent_directory,
            "max_final_t":max_final_t,
            "max_points":max_points, # in time
            "spatial_indices":spatial_indices
            }
        self.initialise_data()
    
    def initialise_data(self):    

        d = self.input_data_args
        target_names,custom_names,sim_output_parent_directory,max_final_t,max_points,spatial_indices = (
            d["target_names"],d["custom_names"],d["sim_output_parent_directory"],d["max_final_t"],d["max_points"],d["spatial_indices"]
        )
        self.target_data = []
        minimum_time_range = np.inf
        lowest_max_points = np.inf           
        for i, mol_name in enumerate(target_names):
            for s_idx in spatial_indices:
                custom_name = custom_names[i]
                dat = PlotData(sim_output_parent_directory,mol_name,"y",max_final_t=max_final_t,max_points=max_points,spatial_index=s_idx,custom_name=custom_name)    
                minimum_time_range = min(dat.get_max_time_range(),minimum_time_range)
                lowest_max_points = min(lowest_max_points, dat.get_num_usable_points())
                self.target_data.append(dat)
        for dat in self.target_data:
            dat.set_max_t(minimum_time_range)
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

    def go(self):
        if not self.check_current():
            self.rerun_ac4dc()
            self.update_outputs()
    
    def initialise_figure(self, plot_title, x_args={}, y_args={}):
        self.fig = go.Figure()
        
        left_anchor = {"xanchor":"left", "x":0.01} 
        #right_anchor = {"xanchor":"right", "x":0.99} 
        right_anchor = {"xanchor":"left", "x":0.6} 
        h_anchor = right_anchor
        title = ""
        if plot_title is not None:
            title = plot_title + " - Free-electron distribution"
        self.fig.update_layout(
            title = title,  # Attention: title overwritten by add_time_slider()
            showlegend= not self.hide_legend,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor = h_anchor["xanchor"],
                x = h_anchor["x"],
                bgcolor = '#F5F5F5',
                font = dict(family="times new roman",size=self.legend_font_size),
                title=self.legend_title,
                title_font=dict(
                    family="times new roman",
                    size=self.font_size,
                ),
            ),               
            # Hide the slider values
            font=dict(
                family="times new roman",
                size=1, 
            )
        )

        self.fig.update_xaxes(x_args)
        self.fig.update_yaxes(y_args)        

    # line_args - a list of kwarg dictionaries to be used for the line argument of go.Scatter().
    # colour_mixup - good for distinguishing plots that are on the same timescale.
    def plot_traces(self, saturation = 0.85, normed = True, colour_mixup = True, line_kwargs = [{},{},{},{},{}], fitE = None):
        Q = 10**(-T_PRECISION)
        # Add a group of traces (each trace corresponding to a point in time) for each target. We will only show one at a time!
        for g, target in enumerate(self.target_data):
            # (Used for colour)
            min_t_c = target.timeData[0]   #
            max_t_c = target.timeData[-1]-(target.timeData[-1]-min_t_c)/len(target.timeData)  # hack, weird indices fix colours and don't affect data.
            # Add traces, one for each slider step
            X = target.energyKnot
            for j, t in enumerate(target.timeData):
                if t//Q*Q > target.max_final_t:
                #if round(t/Q)*Q > target.max_final_t:
                    break

                if j == 0: continue  # Skip t = 0 (empty plot)

                data = target.freeData[j,:]
                if normed:
                    tot = target.get_density(t)
                    data /= tot
                    data/=4*3.14

                ### Cute colours  ★ﾟ~(◠ᴗ ◕✿)ﾟ*｡:ﾟ+ 
                '''                
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
                    t_norm = (t-min_t_c)/(max_t_c-min_t_c)
                    rgb[i] = saturation * rgb_intensity[i] * (1-((t_norm-rgb_bndry[i])/rgb_width[i])**2)
                    rgb[i] = 255*min(1, max(0, rgb[i]))
                    rgba = tuple(rgb) + (a,)
                col = "rgba" + str(tuple(rgba))      
                '''                
                plotly_d3_colors = [
                    '#1f77b4',
                    '#ff7f0e',
                    '#2ca02c',
                    '#d62728',
                    '#9467bd',
                    '#8c564b',
                    '#e377c2',
                    '#7f7f7f',
                    '#bcbd22',
                    '#17becf',
                    ]
                colors = plotly_d3_colors.copy()
                colors.extend(['#4e4e4e'])
                col = colors[g%len(colors)]
                '''
                cmap = plt.cm.get_cmap('viridis') # 'viridis' 'cool' 'plasma' 'inferno' 'cividis'
                if len(self.target_data) == 1:
                    col = cmap(0)
                else:
                    col = cmap(g/(len(self.target_data)-1))
                col = plt.colors.rgb2hex(col)
                '''
                # Choose dependent variable factor depending on if using energy density or electron density.
                density_factor = X # energy density
                if self.use_electron_density:
                    density_factor = 1 
                # Add the trace
                visible = False
                if j == 1:
                    # Make earliest point visible.
                    visible = True
                name=self.target_data[g].target_mol["name"] 
                if self.times_in_legend:
                    name+="%.1f"%t + " fs"
                self.fig.add_trace(
                    go.Scatter(
                        visible=visible,
                        line=dict(color=col, **line_kwargs[g]),
                        name=name,
                        x=X,
                        y=data*density_factor))
                if self.inset:
                    self.fig.add_trace(
                        go.Scatter(
                            visible=visible,
                            line=dict(color=col, **line_kwargs[g]),
                            name=None,
                            x=X,
                            y=data*density_factor,
                            xaxis ='x2',
                            yaxis = 'y2',
                            ))                
                
                    
        # # Plot knots as trace too? (would need to figure out how to make it not move vertically.)
        # knot_to_plot = self.target_data[g].energyKnot
        # self.fig.add_trace(
        #     go.Scatter(
        #         x=knot_to_plot,
        #         y = [10**(self.y_args["range"][0])*1.1]* len(knot_to_plot),
        #         mode="markers",
        #         marker = dict(color='#e66000',size=8),
        #         name="knot"
        #     ),
        # )
        

    #----Widgets----#
    # Time Slider
    def add_time_slider(self,one_slider=False):
        simul_step_slider=False
        if self.num_plots > 1:
            simul_step_slider=True # add a slider that displays all simulations at the same time step. (The limit of plotly: without dash, cannot just turn on or off specific traces via buttons.)
        self.steps_groups = [] # Stores the steps for each plot separately 
        start_step = 0
        time_slider = []
        simul_steps=[]
        for g in range(self.num_plots):
            steps = []
            target = self.target_data[g]
            displaying = "<span style='font-size: 28px; font-family: times new roman'>" +"Displaying:                  </span>"
            subplot_title = dict(text= displaying + "<span style='font-size: "+str(self.font_size)+"px;color:"+ target.title_colour +"; font-family: times new roman'>" + target.target_mol["name"]  + "</span>", yanchor = "top", xanchor = "left", pad = dict(b = 0,l=-400))  # margin-top:100px; display:inline-block;
            allplot_title = copy.deepcopy(subplot_title) 
            allplot_title_colour = "#4d50b3" 
            allplot_title["text"] = displaying + "<span style='font-size: "+str(self.font_size)+"px;color:"+ allplot_title_colour +"; font-family: times new roman'>" + "          All"
            if g == 0:
                self.fig.update_layout({"title": subplot_title})
            for i in range(len(target.timeData) - 1): # -1 as don't have trace for zeroth time step.
                if self.presentation_mode:
                    self.fig.update_layout({"title":"Free electrons"}) 
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
                    trace_label = step["label"]
                    #   Initialise slider that shows all plots.
                    if g == 0:
                        simul_step = copy.deepcopy(step)
                        simul_step["label"] = "  " + trace_label                   
                        simul_steps.append(simul_step)    
                        simul_step["args"][1]["title"] = allplot_title 
                    #   Add later plots' traces at same step (not necessarily same time...).
                    elif self.presentation_mode:
                        # Use time of plot assuming all same.
                        simul_steps[i]["args"][0]["visible"][start_step+i] = True    
                        simul_step["label"] = "  " + trace_label
                        simul_step["args"][1]["title"] = "<span style='font-size: "+str(self.font_size)+"px;color:"+ allplot_title_colour +"; font-family: times new roman'>" + "t = " + trace_label 
                    elif i < len(simul_steps):  
                        simul_steps[i]["args"][0]["visible"][start_step+i] = True    
                        simul_steps[i]["label"] += "  |  " + trace_label
                    if i < len(simul_steps) and g == self.num_plots - 1:
                        # Check if times align. If so, just show one.
                        if simul_steps[i]["label"] == "  " + trace_label  + ("  |  " + trace_label)*g:
                            simul_steps[i]["label"] = trace_label 
            ###
            self.steps_groups.append(steps)
            start_step += len(steps)
        
            # Show individual sliders (Not necessary, can isolate traces by clicking on the legend.)
            individual_sliders = False
            if individual_sliders or len(self.target_data) == 1:
                time_slider.append(dict(
                    active=0,
                    tickwidth=0,
                    tickcolor = "rgba(0,0,0,0)",
                    currentvalue={"prefix": "<span style='font-size: 25px; font-family: times new roman; color = black'>" + target.target_mol["name"] + " - Time [fs]: "},
                    pad={"t": 85+90*g,"r": 200,"l":0},
                    steps=steps,
                    len = 0.5,
                    #font = {"color":"rgba(0.5,0.5,0.5,1)"}
                ))
        if simul_step_slider:
            self.steps_groups.append(simul_steps)
            all_slider = dict(
                active = 0,
                tickwidth = 0,
                currentvalue = {"prefix": "<span style='font-size: 25px; font-family: times new roman; color = white;'>" +"All" + " - Times [fs]: "},
                pad={"t": 85+90*(g+1),"r": 200,"l":0},
                steps = simul_steps,
                #font = {"color":"rgba(0.5,0.5,0.5,1)"}
            )
            if not individual_sliders:
                all_slider = dict(
                    active = 0,
                    tickwidth = 0,
                    currentvalue = {"prefix": "<span style='font-size: "+str(self.font_size)+"px; font-family: times new roman; color = white;'>" +"All" + " - Times [fs]: "},
                    pad={"t": 85,"r": 200,"l":0},
                    steps = simul_steps,
                    len = 0.5,
                    borderwidth = 10,
                    bordercolor = "blue",
                    #font = {"color":"rgba(0.5,0.5,0.5,1)"}
                )                

            if self.presentation_mode:
                all_slider["currentvalue"] =  {"prefix": "<span style='font-size: 100px; font-family: times new roman; color = white;'>" + "t =", "suffix": "<span style='font-size: 100px; font-family: times new roman; color = white;'>" + " fs"}
                #all_slider["pad"] = {"t": 130,"r": 200,"l":0}
                all_slider["pad"] = {"t": -1000,"r": 0,"l":1650}
                all_slider["font"] = {"color":"blue"}
                all_slider["len"] = 0
            time_slider.append(all_slider)
        
        if one_slider:
            time_slider = [time_slider[-1]]
        
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
    def add_scale_button_vertical(self,x_log_args, x_lin_args, y_log_args,y_lin_args,):      
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
                        args=[{'yaxis': y_log_args, 'xaxis': x_lin_args}],
                        label="<br>Log-Lin<br>",
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
                font = {"size": 25,"family": "times new roman"},
        )   
        self.fig.update_layout(
            updatemenus = [scale_button]
        )    
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
                        args=[{'yaxis': y_log_args, 'xaxis': x_lin_args}],
                        label="<br>Log-Lin<br>",
                        method="relayout"    
                    ),    
                    dict(
                        args=[{'yaxis': y_lin_args, 'xaxis': x_lin_args}],
                        label="<br>Lin-Lin<br>",
                        method="relayout"              
                    )
                ]),
                type="buttons",
                # Array buttons along bottom right of graph
                direction = "right",
                pad={"r": 0, "t": 85},
                showactive=True,
                x= 1,
                xanchor="right",
                y= 0, 
                yanchor="top",
                font = {"size": 35,"family": "times new roman"},
        )   
        self.fig.update_layout(
            updatemenus = [scale_button]
        )   

    def get_E_lims(self):
        Emax = -np.inf; Emin = np.inf
        for target in self.target_data:
            Emax = max(Emax,np.max(target.energyKnot))
            Emin = min(Emin,np.min(target.energyKnot))
        return Emin,Emax
            
    # Potential for automatic fitting of MB curve:
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
