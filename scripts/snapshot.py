##### (C) Spencer Passmore 2023 
'''
File: generate_interactive.py
Purpose: For generation of snapshot of electron energy density function.

Usage: Generate snapshots at specified times:
'python3 generate_snapshots lysozyme_3 -15 0 15'

Notes:
By default, assumes that outputs/batches is in AC4DC/output/__Molecular/, and that this file is in AC4DC/scripts

Takes snapshot at time closest to the desired time shared by all simulation output data files. Currently does not throw an error if the specified time is outside the time range 
of one of the sims. TODO make throw error . 


'''
####
import sys
import os
import os.path as path

from types import SimpleNamespace 
from interactive_handler import generate_graphs
from copy import deepcopy

AUTO_VIEW_PLOT = True

HIDE_Y_TICKS = False # in case issue with y ticks causing inconsistent canvas size...


# might have to fiddle a little to get graph canvas to extend to top...

ymax_dict = {
    -12: 0.0012,
    -10: 0.0018,
    0: 0.009,
    10: 0.009,
}

dtick_factor_dict = {
    -12: 2/3,
    -10: 1/2,
    0: 1,
    10: 3*0.010/0.012,
}    

inset_tickformat_dict = {
     -12: ".1f",
    -10: ".1f",
    0: ".0f",
    10: ".0f",         
}


ymax = 0.0024
dtick_factor = 2/3


P = SimpleNamespace()
P.NORMALISE = False
P.ELECTRON_DENSITY = False # if False, use energy densit y  
P.ANIMATION = False # if true, generate animation rather than interactive figure (i.e. automatic slider movement) 
P.ALSO_MAKE_PLOTS = False # Generate static plots for each simulation using _generate_plots.py
P.SINGLE_FRAME = True # Save a png using plotly .
P.NAMING_MODE = 1  # For legend. 0: full details of sim parameters + sim name | 1: elements in sim| 
P.SCALE_DENSITY_BY_THOUSAND = True # Use cubic nm rather than cubic angstrom for measuring energy density  
P.POINTS = 70

width_factor = 0.9
if HIDE_Y_TICKS:
    width_factor *= 0.97
P.SINGLE_FRAME_DICT = dict(
    # In inches
    #width = (7.24409436834*2/3*width_factor),#3.49751*2*3/5,
    width = (7.24409436834*2/3*width_factor),#3.49751*2*3/5,
    #height = 7.24409436834/3, #3.49751*2/3
    height = 7.24409436834/4,#*478/459,
    line_width = 2,
    font_size = 10,
    superscript_font_size = 8,
    
    pad=100,

    show_legend = False, 

    ylog = False, # < uses linear scale if False 
    xlog = False, # <

    #times = [-10],#[-10],#[0],#[-10,0], # fs.  Use multiple to plot traces from multiple times (colours will be same between frames however.)
    #y_range = [0,0.0014],#lin:[0,0.0014]/[0,0.009],#log: [-5,-1],      # < If axis is log scale, limits correspond to powers of 10.
    
    #times = [-10],
    #y_range = [0,0.0014], # -10 fs
    #y_range = [0,0.0025], # -10 fs solvated
    #y_range = [0,0.001], # -10 fs intermediate
    #times = [0],
    #y_range = [0,0.009], # 0 fs
    #y_range = [0,0.014], # 0 fs solvated
    y_range = [0,ymax],
    x_range = [None,7500], # (use [None,None] for default)
)
P.INSET = True
P.SUBPLOT_AX_FONT_COLOUR = "#495187"# "blue"
P.SUBPLOT_AX_GRID_COLOUR = "#c4c4eb"
P.INSET_DICT = dict(
    axes_kwargs = dict(
        xaxis2 =dict(
            #domain=[0.1, 0.4],
            domain=[0.07, 0.37],
            range = [0,120], #[0,50],
            anchor='y2',
            gridcolor=P.SUBPLOT_AX_GRID_COLOUR, 
            zeroline = False,
            #zerolinecolor = P.SUBPLOT_AX_GRID_COLOUR
        ),
        yaxis2=dict(
            #domain=[0.15, 0.8],
            domain=[0.11, 0.76],
            range = [0,None],
            #range = [0,2.5*P.SINGLE_FRAME_DICT['y_range'][1]],
            anchor='x2',
            gridcolor=P.SUBPLOT_AX_GRID_COLOUR, 
            zeroline = False,
            tickformat = ".1f",
            tick0=P.SINGLE_FRAME_DICT['y_range'][0],
            dtick = (P.SINGLE_FRAME_DICT['y_range'][1]- P.SINGLE_FRAME_DICT['y_range'][0])/3*2*(1+999*P.SCALE_DENSITY_BY_THOUSAND),          
            #dtick = (P.SINGLE_FRAME_DICT['y_range'][1]- P.SINGLE_FRAME_DICT['y_range'][0])/3*3*(1+999*P.SCALE_DENSITY_BY_THOUSAND),          
            #dtick = (P.SINGLE_FRAME_DICT['y_range'][1]- P.SINGLE_FRAME_DICT['y_range'][0])/3*6*(1+999*P.SCALE_DENSITY_BY_THOUSAND),          
            #zerolinecolor = P.SUBPLOT_AX_GRID_COLOUR
        ),
    ),
)
if HIDE_Y_TICKS:
    P.INSET_DICT['axes_kwargs']['yaxis']=dict(showticklabels=False)
    #P.SINGLE_FRAME_DICT['axes_kwargs'] = {}
    #P.SINGLE_FRAME_DICT['axes_kwargs']['yaxis']=dict(showticklabels=False)

if  len(sys.argv) < 3:
    print("Usage: Generate snapshots at some number of times, e.g. t = -10 and t = 0, with: 'python3 scripts/"+path.basename(__file__)+" lysozyme_3 tetrapeptide_1 -10 0 10'")
    exit()
n=0
for k in sys.argv:
   if  k.strip('-').isnumeric() or len(k.strip('-').split('.'))==2 and k.strip('-').split('.')[0].join(k.strip('-').split('.')[1]).isnumeric() :
       break 
   n+=1
   
def get_latest_plot_name():
    dir = path.abspath(path.join(__file__ ,"../../output/_Graphs/plots")) + "/" 
    import glob
    files = glob.glob(dir+"*")
    latest_file = max(files,key=path.getctime)
    return latest_file
    
def view_plot(plot_name):
    os.system("wslview " + plot_name)
    
def update_dict(t):
    P.SINGLE_FRAME_DICT["times"] = [t]  # Put at value to cutoff times early.
    P.SINGLE_FRAME_DICT['y_range'] = [0,ymax_dict[t]]
    P.INSET_DICT['axes_kwargs']['yaxis2']['dtick'] = dtick_factor_dict[t]*(P.SINGLE_FRAME_DICT['y_range'][1]- P.SINGLE_FRAME_DICT['y_range'][0])*(1+999*P.SCALE_DENSITY_BY_THOUSAND)
    P.INSET_DICT['axes_kwargs']['yaxis2']['tickformat'] = inset_tickformat_dict[t]
    
if n >= len(sys.argv):
    print("Times not provided.")
    exit()
plotted_files = []
file_name_tag=""
if not HIDE_Y_TICKS:
    file_name_tag="withAxis"   # TODO temporary for my convenience
for snapshot_t in sys.argv[n:]:
    print("Taking snapshot at t = "+snapshot_t+" fs")
    #P.END_T = float(snapshot_t)  # Put at value to cutoff times early.
    update_dict(float(snapshot_t))
    generate_graphs(deepcopy(P),sys_argv = sys.argv[0:n],file_name_tag=file_name_tag) # deepcopy because plotly modifies the dicts
    plotted_files.append(get_latest_plot_name())
if AUTO_VIEW_PLOT:
    for plot in plotted_files:
        view_plot(plot)

############### Scratchpad
# python3.9 scripts/generate_snapshot.py lys_nass_no_S_3 lys_nass_gauss lys_nass_Gd_gauss_1 


# python3.9 scripts/snapshot.py  lys_solvated_light_4 lys_solvated_6 lys_Gd_solvated_1 lys_Gd_salt_solvated_1 -10 0 10

# lys_solvated_light_H_9 lys_solvated_H_2 Gd_no_salt_2 Gd_salt_5
# no_S_1 no_Gd_and_no_salt_3 Gd_and_no_salt_3 Gd_and_salt_1

# python3.9 scripts/snapshot.py lys_solvated_light_10  lys_solvated_16 CB_no_salt_6 CB_11 -10 0 10
