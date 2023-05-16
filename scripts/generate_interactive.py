##### (C) Spencer Passmore 2023 
'''
File: generate_interactive.py
Purpose: For generation of interactive electron density figure, alongside misc plots.

1. generate figures for single output
'python3 generate_interactive.py lysozyme_3'

2. Compare interactive figures 
'python3 generate_interactive.py lysozyme_3 lysozyme_4'
If two handles provided dotted line is used for first.  

3. Use a folder containing a "batch" of simulation output folders to generate a folder of graphs (batch folder in /output/__Molecular/)
'python3 generate_interactive.py batch_lys'  
if made with generate_batch.py, batch_lys likely contains files like: lys-1_2, lys-2_2, lys-3_2,... but names do not matter.

TODO ac4dc doesnt auto put results into a batch folder currently.
TODO need to be able to specify folders containing outputs e.g. to compare a few outputs
TODO add -N flag for normalisation (cld also add normalisation option to plotly to replace sliders with normed trace vrsns.)
Notes:
Simulation outputs/batches should be in AC4DC/output/__Molecular/ while this scripts is in AC4DC/scripts


'''
####


terminal_mode = True
normalise = False
xmin, xmax = 1, 2e4

END_T = 99#10 #0.0
POINTS = 70
#####

import os
import os.path as path
from interactive_core import InteractivePlotter#, fit_maxwell, maxwell
from _generate_plots import make_some_plots
import sys
import numpy as np
from QoL import set_highlighted_excepthook
sys.path.append(os.path.join(os.path.dirname(__file__), 'scattering'))
from core_functions import get_sim_params # so that we can name batches for the parameters we care about. 

def main():
    set_highlighted_excepthook()
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    graph_folder = path.abspath(path.join(__file__ ,"../../output/_Graphs/")) + "/" 
    sim_input_dir = path.abspath(path.join(__file__,"../../input/")) + "/"

    if terminal_mode:
        sys_argv = sys.argv  
    if  len(sys_argv) < 2:
        print("Usage: python3 " +path.basename(__file__)+ " Carbon_1 <comparison_handle1...etc. (optional)>. Alternatively a folder name can be provided to generate interactives for all contained outputs.")
        exit()
   
    def make_outfolder(outdir,subfolders):
        for subdir in subfolders:
            try:
                os.makedirs(outdir+subdir,exist_ok=True)
            except OSError:
                print("Ignoring os.makedirs error.")
                pass
    
    def single_interactive(target_handles,data_parent_dir,outdir):
        print(target_handles,data_parent_dir,outdir)
        custom_names = get_custom_names(target_handles,data_parent_dir)
        fname_out = target_handles[0]
        if len(sys_argv) > 2:
            fname_out += '_' + sys_argv[2][:15] + "..-"
        plot_title = fname_out.replace('_',' ')  
        int_subdir = "interactives/"; plot_subdir = "plots/"  # Separate out types of plots
        make_outfolder(outdir,[int_subdir,plot_subdir])
        generate(target_handles,data_parent_dir,plot_title,fname_out + "_Int", normalise,outdir+int_subdir,custom_names)    
        for i, target in enumerate(target_handles):
            if len(sys_argv) > 2:
                fname_out = target
            make_some_plots(target,data_parent_dir,fname_out+"_Plt",outdir+plot_subdir,True,False,True,False)   #TODO turning off bound plots for now, they are joined with multiple species.          
    def get_custom_names(target_handles,sim_data_parent_dir):
        custom_names = []
        for handle in target_handles:
            sim_params = list(get_sim_params(sim_input_dir,sim_data_parent_dir,handle))
            sim_params = sim_params[2:] # 0 and 1 are start and end time,
            unit_dict = sim_params.pop()
            param_dict = sim_params.pop()
            name = ["  " + param_dict[i][0] + ": " +str(p) + unit_dict[i]  for i, p in enumerate(sim_params)]
            name = ''.join(name) + " [" + handle + "]"
            name = name[2:]  
            custom_names.append(name)
        return custom_names

    given_path = ""
    is_parent_dir = False
    if len(sys_argv) == 2:
        is_parent_dir = True
        given_path = molecular_path + sys_argv[1] + "/"
        if path.isfile(given_path):
            raise Exception("Please pass simulation output folder name or folder containing simulation output folders")
        for f in os.listdir(given_path):
            if path.isfile(given_path + "/"+ f):
                is_parent_dir = False
                break
    if is_parent_dir:
        # Generate interactive for each output contained in the parent directory
        outdir = graph_folder + sys_argv[1] + "/"
        for target_handle in os.listdir(given_path):
            single_interactive([target_handle],given_path,outdir)
        # ipl = multi_interactive(given_path)
        #save_interactive(ipl,outdir,sys_argv[1] + "_interactive")
        print("Interactives successfully saved to",outdir)
    else: 
        outdir = graph_folder
        target_handles = []  
        # Store folder names from commandline args without .mol extension. 
        for i, target in enumerate(sys_argv):
            if i == 0: continue #####
            if target.endswith(".mol"):
                target = target[0:-4]         
            target_handles.append(target)        
        single_interactive(target_handles,molecular_path,outdir)
    print("Done! Remember to drink water!")

# Not used currently, from https://gist.github.com/rayliuca/a50f86acaae019b16e1621e81c2cf4ea
def symlog(arr, base=10, linthresh=1e-10, linscale=1):
    _linscale_adj = (linscale / (1.0 - base ** -1))
    _log_base = np.log(base)
    arr = np.array(arr)
    abs_arr = np.abs(arr)
    
    # arrays with all elements within the linear region of this symlog transformation 
    linear = np.max(abs_arr, axis=0) < linthresh

    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.sign(arr) * linthresh * (
            _linscale_adj +
            np.log(abs_arr / linthresh) / _log_base)
        inside = abs_arr <= linthresh

    out[inside] = arr[inside] * _linscale_adj

    out = out/np.max(np.abs(out), axis=0)*np.log(np.max(abs_arr, axis=0))/_log_base

    out[np.array([linear]*out.shape[0])] = arr[np.array([linear]*out.shape[0])]
    
    return out



def set_up_interactive_axes(ipl,plot_title):
    # Y scaling
    lin_ymax = 0.03
    lin_ymin = -0.005
    log_ymax = 1
    log_ymin = 1e-4      
    # Axis parameters, these define axis properties depending on what scale is used.
    xlabel = 'Energy (eV)'
    ylabel = '$f(\\epsilon) \\Delta \\epsilon$'   #TODO: Get this working on offline file saves somehow.
    x_log_args = {'title': {"text": xlabel + " - log scale", "font":{"size": 30,"family": "roboto"}}, 'tickfont': {"size": 20}, 'type' : "log", "range" : [np.log10(xmin),np.log10(xmax)]}
    x_lin_args = {'title': {"text": xlabel + " - lin scale", "font":{"size": 30,"family": "roboto"}}, 'tickfont': {"size": 20}, 'type' : "linear", "range" : [xmin,xmax]}
    y_lin_args = {'title': {"text": ylabel + " - lin scale", "font":{"size": 30,"family": "roboto"}}, 'tickfont': {"size": 20}, 'type' : "linear", "range" : [lin_ymin,lin_ymax]}
    y_log_args = {'title': {"text": ylabel + " - log scale", "font":{"size": 30,"family": "roboto"}}, 'tickfont': {"size": 20}, 'type' : "log", "range" : [np.log10(log_ymin),np.log10(log_ymax)]}
    
    # Initialises graph object.
    ipl.initialise_interactive(plot_title, x_log_args,y_log_args) 
    return (x_log_args,x_lin_args,y_log_args,y_lin_args)

def generate(target_handles,sim_data_parent_dir,plot_title,fname_out,normalise,outdir,custom_names = None):
    '''
    Arguments:
    target_handles: The name of the folder containing the simulation's data (the csv files). 
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''
    print("[ Interactive ] Creating electron density interactive figure")
    # Initialises plotter object with data from files.
    ipl = InteractivePlotter(target_handles,sim_data_parent_dir, max_final_t=END_T,max_points=POINTS,custom_names=custom_names)  
    scale_button_args = set_up_interactive_axes(ipl,plot_title)
    # The meat of the plotting. Plot line for each point in time
    line_1 = {'width': 6,"dash": '10,1'}
    line_2 = {'width': 6,}
    line_kwargs = [line_1,line_2]
    if len(target_handles) != 2:
        line_kwargs = [line_2 for l in range(len(target_handles))]
    ipl.plot_traces(normed=normalise, line_kwargs=line_kwargs)
    # Add widgets
    ipl.add_scale_button(*scale_button_args)                 
    ipl.add_time_slider()  
    #-----Save-----#
    save_interactive(ipl,outdir,fname_out)

def save_interactive(ipl,outdir,fname_out):
    extension = ".html"
    file_path =  outdir + fname_out + extension
    ipl.fig.write_html(file_path)

    print("Done!")

if __name__ == "__main__":
    main()



# <Would need to use dash implementation for this:>.
# def multi_interactive(output_parent_directory):
#     target_handles = []
#     for target_handle in os.listdir(given_path):
#         target_handles.append(target_handle)
#     ipl = InteractivePlotter(target_handles,output_parent_directory, max_final_t=END_T,max_points=POINTS) 
#     add_params_to_multi_interactive(ipl,target_handles)
#     return  

# def add_params_to_multi_interactive(ipl,target_handles):
#     for idx,target_handle in enumerate(target_handles):
#         energy, fwhm, photons = [int(val) for val in target_handle.split("_")[:3]]
#         ipl.set_trace_params(idx,energy,fwhm,photons,normed=normalise)
#     ipl.plot_traces()
#     ipl.switching_time_sliders()
