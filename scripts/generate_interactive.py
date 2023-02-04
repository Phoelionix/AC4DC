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

    set_highlighted_excepthook()

    if terminal_mode:
        sys_argv = sys.argv  
    if  len(sys_argv) < 2:
        print("Usage: python3 _save_interactive.py Carbon_1 <name_suffix (optional)>")
        exit()        

    target_folder_names = []  
    # Store folder names from commandline args without .mol extension. 
    for i, target in enumerate(sys_argv):
        if i == 0: continue #####
        if target.endswith(".mol"):
            target = target[0:-4]         
        target_folder_names.append(target)

    # Naming file/plot
    fname_out = target_folder_names[0]
    if len(sys_argv) > 2:
        fname_out += '_' + sys_argv[2][:10] + "..-"
    fname_out += "_interactive"
    plot_title = fname_out.replace('_',' ')
    
    # Y scaling
    lin_ymax = 0.035
    lin_ymin = -0.035
    log_ymax = 1
    log_ymin = 1e-10      
    normalise = True  # TODO make True default, and False an option for cmdline arg.
    # Axis params
    xlabel = 'Energy (eV)'
    ylabel = '$f(\\epsilon) \\Delta \\epsilon$'   #TODO: Get this working on offline file saves somehow.
    x_log_args = {'title': xlabel + " - log scale", 'type' : "log", "range" : [np.log10(xmin),np.log10(xmax)]}
    x_lin_args = {'title': xlabel + " - lin scale", 'type' : "linear", "range" : [xmin,xmax]}
    y_lin_args = {'title': ylabel + " - lin scale", 'type' : "linear", "range" : [lin_ymin,lin_ymax]}
    y_log_args = {'title': ylabel + " - log scale", 'type' : "log", "range" : [np.log10(log_ymin),np.log10(log_ymax)]}
    #
    # Initialises plotter object with data from files.
    ipl = InteractivePlotter(target_folder_names,"y")  
    # Initialises graph object.
    ipl.initialise_interactive(plot_title, x_log_args,y_log_args) 
    # The meat of the plotting. Plot line for each point in time
    ipl.plot_traces(normed=normalise)
    # Add widgets
    ipl.add_scale_button(x_log_args,x_lin_args, y_log_args,y_lin_args)                 
    ipl.add_time_slider()  
    #ipl.add_simulation_menu()

    # from dash import Dash, html, Input, Output
    # from dash import dcc
    # import dash_daq as daq

    # app = Dash(__name__)

    # # app.layout = html.Div([
    # #     daq.Slider(
    # #         id='my-daq-slider-ex-1',
    # #         value=17
    # #     ),
    # #     html.Div(id='slider-output-1')
    # # ])

    # app.layout = html.Div([
    #     dcc.Graph(
    #         figure=ipl.fig,
    #         id='my-daq-slider-ex-1',
    #     ),
    #     html.Div(id='slider-output-1')
    # ])

    # @app.callback(
    #     Output('slider-output-1', 'children'),
    #     Input('my-daq-slider-ex-1', 'value')
    # )
    # def update_output(value):
    #     return f'The slider is currently at {value}.'
    # app.run_server(debug=True)
    # #ipl.plot_maxwell(44.1,0.06*3/2)  # ~Sanders -7.5 fs - density of MB assumed to be 50% of total.  
    # #ipl.fig.show()
    #-----Save-----#
    extension = ".html"
    outdir = "../../../AC4DC_Interactives/"
    file_path = path.abspath(path.join(__file__ ,outdir + fname_out + extension))
    ipl.fig.write_html(file_path)

    print("Done!")

if __name__ == "__main__":
    main()
# %%
