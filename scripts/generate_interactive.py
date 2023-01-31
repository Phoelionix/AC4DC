#####
terminal_mode = True
xmin, xmax = 1, 1e4
#####

import os.path as path
from interactive_core import InteractivePlotter#, fit_maxwell, maxwell
import sys
import numpy as np
from QoL import set_highlighted_excepthook

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

def main():
    print("[ Interactive ] Creating electron density interactive figure")

    if terminal_mode:
        sys_argv = [sys.argv[0],sys.argv[1]]    
    target_folder_name = sys_argv[1]     
    if target_folder_name.endswith(".mol"):
        target_folder_name = target_folder_name[0:-4] 

    set_highlighted_excepthook()

    # Basic num arguments check
    if  len(sys_argv) < 2:
        print("Usage: python3 _save_interactive.py Carbon_1 <name_suffix (optional)>")
        exit()

    # Naming file/plot
    label = target_folder_name
    if len(sys_argv) > 2:
        label += '_' + sys_argv[2]
    label += "_interactive"
    name = target_folder_name.replace('_',' ')
    

    # Initialisation
    lin_ymax = 0.035
    lin_ymin = -0.035
    log_ymax = 1
    log_ymin = 1e-10      
    normalise = True
    # Axis params
    xlabel = 'Energy (eV)'
    ylabel = '$f(\\epsilon) \\Delta \\epsilon$'   #TODO: Get this working on offline file saves somehow.
    x_args = {'title': xlabel, 'type' : "log", "range" : [np.log10(xmin),np.log10(xmax)]}
    y_lin_args = {'title': ylabel + " (lin)", 'type' : "linear", "range" : [lin_ymin,lin_ymax]}
    y_log_args = {'title': ylabel + " (log)", 'type' : "log", "range" : [np.log10(log_ymin),np.log10(log_ymax)]}
    #
    ipl = InteractivePlotter(target_folder_name,"y")
    ipl.initialise_interactive(name, x_args,y_log_args)
    # Plot line for each point in time
    ipl.plot_traces(normed=normalise)
    # Add widgets
    ipl.add_scale_button(y_log_args,y_lin_args)                 
    ipl.add_time_slider()  

    #ipl.plot_maxwell(44.1,0.06*3/2)  # ~Sanders -7.5 fs - density of MB assumed to be 50% of total.  
    #ipl.fig.show()
    #-----Save-----#
    extension = ".html"
    outdir = "../../../AC4DC_Interactives/"
    file_path = path.abspath(path.join(__file__ ,outdir + label + extension))
    ipl.fig.write_html(file_path)

    print("Done!")

if __name__ == "__main__":
    main()
# %%
