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
from core_functions import PlotData, get_mol_file, parse_elecs_from_latex, get_spatial_indices, ATOMS, ATOMNO, T_PRECISION
import matplotlib as plt
import inspect
pio.templates.default = "seaborn" #"plotly_dark" # "plotly"


class SpatialInput:
   # max_final_t, float, end time in femtoseconds. Not equivalent to time duration
    # max_points, int, number of points (within the timespan) for the interactive to have at maximum.
    def __init__(self, target_names, sim_output_parent_directory, spatial_indices, max_final_t, max_points, custom_names, use_electron_density):
        '''
        output_parent_directory: absolute path
        use_electron_density: If True, plot electron density rather than energy density
        Various changes, such as bigger font, etc. for presentation purposes.
        '''
        

        max_points+=1  # (we exclude t=0)
        self.multi_trace_params = [""]*len(target_names)  # 
        self.use_electron_density = use_electron_density
        self.presentation_mode = False

        # self.times_in_legend = times_in_legend
        # self.hide_legend = hide_legend
        # self.legend_title = legend_title

        if spatial_indices is None:
            spatial_indices = get_spatial_indices( sim_output_parent_directory, target_names[0])
            for target in target_names:
                assert(spatial_indices==get_spatial_indices( sim_output_parent_directory, target))            
        
        self.num_cells = len(spatial_indices)
        self.num_plots = len(spatial_indices)*len(target_names)
        if custom_names is None:
            custom_names = [None]*self.num_plots
        self.data_args = {
            "target_names":target_names,
            "custom_names":custom_names,
            "sim_output_parent_directory":sim_output_parent_directory,
            "max_final_t":max_final_t,
            "max_points":max_points,
            "spatial_indices":spatial_indices
            }
        
class SpatialData:
    def __init__(self, input : SpatialInput ):
        self.input_args = input
        self.initialise_data(input.data_args)

    def initialise_data(self, d : dict):    

        target_names,custom_names,sim_output_parent_directory,max_final_t,max_points,spatial_indices = (
            d["target_names"],d["custom_names"],d["sim_output_parent_directory"],d["max_final_t"],d["max_points"],d["spatial_indices"]
        )    
        target_data = []
        minimum_time_range = np.inf
        lowest_max_points = np.inf           
        for i, mol_name in enumerate(target_names):
            for s_idx in spatial_indices:
                custom_name = custom_names[i]
                dat = PlotData(sim_output_parent_directory,mol_name,"y",max_final_t=max_final_t,max_points=max_points,spatial_index=s_idx,custom_name=custom_name)    
                minimum_time_range = min(dat.get_max_time_range(),minimum_time_range)
                lowest_max_points = min(lowest_max_points, dat.get_num_usable_points())
                target_data.append(dat)
        for dat in target_data:
            dat.set_max_t(minimum_time_range)
            dat.max_points = lowest_max_points
            dat.update_outputs()        
        self.target_data = target_data

    
class SpatialInteractive:
    def __init__(self, target_names, sim_output_parent_directory, spatial_indices = None, max_final_t = 30, max_points = 70, custom_names = None,use_electron_density = False):
        tmp_dict = inspect.getargvalues(inspect.currentframe())[-1]
        tmp_dict.pop('self')
        print(tmp_dict)
        self.input = SpatialInput(**tmp_dict)
        self.target_data = SpatialData(self.input).target_data
        self.new_fig()

    def new_fig(self,width=800,height=800):
        self.fig = go.Figure()

        # Set figure size
        self.fig.update_layout(width=width, height=height)   

    def save_interactive(self,outdir,fname_out):
        extension = ".html"
        file_path =  outdir + fname_out + extension
        self.fig.write_html(file_path)

        print("Done!") 

    def plot_circles(self, atom : str, thickness : float, R : float):
        min_z = 0
        max_z = 1

        self.fig.update_xaxes(range=[-R*1.1, R*1.1], zeroline=False)
        self.fig.update_yaxes(range=[-R*1.1, R*1.1])

        # Build up 'layered cake' of circles starting from base
        assert self.input.num_cells==self.input.num_plots  # for now
        for i in reversed(range(len(self.target_data))):
            r = (1+i)*thickness
            target = self.target_data[i]
            assert target.spatial_index == i
            charge = target.get_charge(atom)
            print(f"Average charge:{charge[-1]}")
            depletion=charge/ATOMNO[atom]

            # colours
            col_z = depletion
            col_z -= min_z
            col_z /=(max_z-min_z)
            cmap = plt.cm.get_cmap('magma')



            Q = 10**(-T_PRECISION)
            # Add traces, one for each slider step
            for j, t in enumerate(target.timeData):
                if t//Q*Q > target.max_final_t:
                #if round(t/Q)*Q > target.max_final_t:
                    break

                # Start with earliest point in time visible (this is overridden by call to time_slider).
                #visible = True
                # if j == 0:
                #     visible = True
                name=target.target_mol["name"]+str(target.spatial_index) 
                # if self.times_in_legend:
                #     name+="%.1f"%t + " fs"

                col = cmap(col_z[j])
                col = plt.colors.rgb2hex(col)


                self.fig.add_shape(type="circle",
                    visible=True,
                    opacity=1,
                    xref="x", yref="y",
                    fillcolor=col,#"PaleTurquoise",
                    x0=-r, y0=-r, x1=r, y1=r,
                    line_color=col,#"LightSeaGreen",
                    name=name,
                )
        self.add_time_slider(True)
        self.fig.update_layout(yaxis_scaleanchor="x")


    #----Widgets----#
    # Time Slider
    def add_time_slider(self,one_slider=False,font_size=30):

        simul_step_slider=False
        if self.input.num_plots > 1:
            simul_step_slider=True # add a slider that displays all simulations at the same time step. (The limit of plotly: without dash, cannot just turn on or off specific traces via buttons.)
        self.steps_groups = [] # Stores the steps for each plot separately 
        start_step = 0
        time_slider = []
        simul_steps=[]
        for g in range(self.input.num_plots):
            steps = []
            target = self.target_data[g]
            displaying = "<span style='font-size: 28px; font-family: times new roman'>" +"Displaying:                  </span>"
            subplot_title = dict(text= displaying + "<span style='font-size: "+str(font_size)+"px;color:"+ target.title_colour +"; font-family: times new roman'>" + target.target_mol["name"]  + "</span>", yanchor = "top", xanchor = "left", pad = dict(b = 0,l=-400))  # margin-top:100px; display:inline-block;
            allplot_title = copy.deepcopy(subplot_title) 
            allplot_title_colour = "#4d50b3" 
            allplot_title["text"] = displaying + "<span style='font-size: "+str(font_size)+"px;color:"+ allplot_title_colour +"; font-family: times new roman'>" + "          All"
            if g == 0:
                self.fig.update_layout({"title": subplot_title})
            for i in range(len(target.timeData)):
                if self.input.presentation_mode:
                    self.fig.update_layout({"title":"Free electrons"}) 
                
                if  i <len(target.timeData)-1 and target.timeData[i+1] > target.max_final_t:
                    break                
                step = dict(
                    method="relayout",
                    #method="update",
                    args=[
                        #{"visible": [False] * len(self.fig.data)},  # style attribute
                        {"title": subplot_title, "shapes": [self.fig.layout.shapes[start_step + i]]},                   # layout attribute  
                        #{"title": subplot_title,},                   # layout attribute         #Add line below: + '<br>' +  '<span style="font-size: 12px;">line2</span>'}  
                    ],      
                    label= "  " + "%.2f" % target.timeData[i],
                )
                #step["args"][0]["shapes"][start_step + i]["visible"] = True  #  When at this step toggle i'th shape in target's group to "visible"
                #step["args"][0]["visible"][start_step + i] = True  #  When at this step toggle i'th trace in target's group to "visible"
                steps.append(step)
                if simul_step_slider:
                    trace_label = step["label"]
                    #   Initialise slider that shows all plots.
                    if g == 0:
                        simul_step = copy.deepcopy(step)
                        simul_step["label"] = "  " + trace_label                   
                        simul_steps.append(simul_step)    
                        simul_step["args"][0]["title"] = allplot_title 
                    #   Add later plots' traces at same step (not necessarily same time...).
                    elif self.input.presentation_mode:
                        # Use time of plot assuming all same.
                        #simul_steps[i]["args"][0]["visible"][start_step+i] = True    
                        simul_steps[i]["args"][0]["shapes"].extend(steps[i]["args"][0]["shapes"])
                        simul_step["label"] = "  " + trace_label
                        simul_step["args"][0]["title"] = "<span style='font-size: "+str(font_size)+"px;color:"+ allplot_title_colour +"; font-family: times new roman'>" + "t = " + trace_label 
                    elif i < len(simul_steps):  
                        #simul_steps[i]["args"][0]["visible"][start_step+i] = True    
                        simul_steps[i]["args"][0]["shapes"].extend(steps[i]["args"][0]["shapes"])
                        simul_steps[i]["label"] += "  |  " + trace_label
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
                    currentvalue = {"prefix": "<span style='font-size: "+str(font_size)+"px; font-family: times new roman; color = white;'>" +"All" + " - Times [fs]: "},
                    pad={"t": 85,"r": 200,"l":0},
                    steps = simul_steps,
                    len = 0.5,
                    borderwidth = 10,
                    bordercolor = "blue",
                    #font = {"color":"rgba(0.5,0.5,0.5,1)"}
                )                

            if self.input.presentation_mode:
                all_slider["currentvalue"] =  {"prefix": "<span style='font-size: 100px; font-family: times new roman; color = white;'>" + "t =", "suffix": "<span style='font-size: 100px; font-family: times new roman; color = white;'>" + " fs"}
                #all_slider["pad"] = {"t": 130,"r": 200,"l":0}
                all_slider["pad"] = {"t": -1000,"r": 0,"l":1650}
                all_slider["font"] = {"color":"blue"}
                all_slider["len"] = 0
            time_slider.append(all_slider)
        
        if one_slider: # if multiple plots, this changes the time slider to the simulstep_slider
            time_slider = [time_slider[-1]]
        

        self.fig.update_layout(sliders=time_slider)   
        

        self.fig.update_layout(**time_slider[-1]["steps"][0]["args"][0])   # set layout to earliest time point
        
