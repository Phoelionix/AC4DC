##### (C) Spencer Passmore 2023 
terminal_mode = True
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

def main():
    set_highlighted_excepthook()
    normalise = True
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"

    if terminal_mode:
        sys_argv = sys.argv  
    if  len(sys_argv) < 2:
        print("Usage: python3 " +path.basename(__file__)+ " Carbon_1 <comparison_handle1...etc. (optional)>. Alternatively a folder name can be provided to generate interactives for all contained outputs.")
        exit()        
   
    def make_outfolder(outdir):
        try:
            os.makedirs(outdir)
        except OSError:
            print("Ignoring os.makedirs error.")
            pass
    
    def single_interactive(target_handles,output_parent_directory,outdir):
        fname_out = target_handles[0]
        if len(sys_argv) > 2:
            fname_out += '_' + sys_argv[2][:15] + "..-"
        plot_title = fname_out.replace('_',' ')   
        make_outfolder(outdir)
        generate(target_handles,plot_title,fname_out + "_interactive", normalise,output_parent_directory,outdir)    
        for target in target_handles:
            make_some_plots(target,output_parent_directory,fname_out+"_plots",outdir,True,True,True,False)            
    
    # <Would need to use dash implementation>.
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

    given_path = ""
    is_parent_dir = True
    if len(sys_argv) < 3:
        given_path = molecular_path + sys_argv[1] + "/"
        for f in os.listdir(given_path):
            if path.isfile(given_path + "/"+ f):
                is_parent_dir = False
                break
        # Generate interactive for each output contained in the parent directory
    if is_parent_dir:
        outdir = path.abspath(path.join(__file__ ,"../../../AC4DC_Interactives/")) + "/" + sys_argv[1] + "/"
        for target_handle in os.listdir(given_path):
            single_interactive([target_handle],given_path,outdir)
        # ipl = multi_interactive(given_path)  
        #fname_out = sys_argv[1] + "_interactive"
        #save_interactive(ipl,outdir,fname_out)
        print("Interactives successfully saved to",outdir)
    else: 
        outdir = path.abspath(path.join(__file__ ,"../../../AC4DC_Interactives/")) + "/"
        target_handles = []  
        # Store folder names from commandline args without .mol extension. 
        for i, target in enumerate(sys_argv):
            if i == 0: continue #####
            if target.endswith(".mol"):
                target = target[0:-4]         
            target_handles.append(target)        
        print("ADAS",target_handles)
        single_interactive(target_handles,molecular_path,outdir)


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

def generate(target_handles,plot_title,fname_out,normalise,output_parent_directory,outdir):
    '''
    Arguments:
    target_handles: The output folder name (not directory)
    output_parent_directory: absolute path
    '''

    print("[ Interactive ] Creating electron density interactive figure")
    # Initialises plotter object with data from files.
    ipl = InteractivePlotter(target_handles,output_parent_directory, max_final_t=END_T,max_points=POINTS)  
    scale_button_args = set_up_interactive_axes(ipl,plot_title)
    # The meat of the plotting. Plot line for each point in time
    line_1 = {'width': 6,"dash": '10,1'}
    line_2 = {'width': 6,}
    line_kwargs = [line_1,line_2]
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
# %%
