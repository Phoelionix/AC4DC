#%%
import numpy as np
from script_frame_plot_core import Plotter
from IPython.display import clear_output
import plotly.io as pio
import plotly.offline as plot_off
import IPython
import time

def plot_frame(ipl,energies,densities,knot_to_plot):

    ipl.update_data(energies,densities)
    normalise = True # Currently we use fixed axes.
    # Plot the trace of the electron energy density (f[e]*e) 
    ipl.plot_step(normed=normalise)
    # Plot the knots (the grid/basis points)
    ipl.plot_the_knot(knot_to_plot)   


def plot_frame_from_c(energies,densities,knot_points):
    knot_points = [27.211385*e for e in knot_points]  # knots are passed in Ha, sorry.
    pio.renderers.default = "notebook"
    # Axis params
    log_ymin, log_ymax  = 1e-4, 1      
    xmin, xmax = 1, energies[np.max(np.nonzero(energies)) - 1]   # from 1 eV to max ignoring upper boundary knot.
    xlabel = 'Energy (eV)'
    ylabel = '$f(\\epsilon) \\Delta \\epsilon$'   #TODO: Get this working on offline file saves somehow.
    x_log_args = {'title': {"text": xlabel + " - log scale", "font":{"size": 30,"family": "roboto"}}, 'tickfont': {"size": 20}, 'type' : "log", "range" : [np.log10(xmin),np.log10(xmax)]}
    y_log_args = {'title': {"text": ylabel + " - log scale", "font":{"size": 30,"family": "roboto"}}, 'tickfont': {"size": 20}, 'type' : "log", "range" : [np.log10(log_ymin),np.log10(log_ymax)]}
    #    
    ipl = Plotter()
    ipl.initialise_figure("Current Distribution", x_log_args,y_log_args) 
    plot_frame(ipl,energies,densities,knot_points)
    #plot_off.iplot(ipl.fig)
    #ipl.fig.show()
    #IPython.display.display(ipl.fig);     
    pio.write_image(ipl.fig, "_live_plot.png");  
def clear_frame_c():
    clear_output(wait=True)

if __name__ == "__main__":
    # Test whether frames are shown properly
    import time
    # Axis params
    log_ymin, log_ymax  = 1e-4, 1      
    xmin, xmax = 1, 1e4    
    xlabel = 'Energy (eV)'
    ylabel = '$f(\\epsilon) \\Delta \\epsilon$'   #TODO: Get this working on offline file saves somehow.
    x_log_args = {'title': {"text": xlabel + " - log scale", "font":{"size": 30,"family": "roboto"}}, 'tickfont': {"size": 20}, 'type' : "log", "range" : [np.log10(xmin),np.log10(xmax)]}
    y_log_args = {'title': {"text": ylabel + " - log scale", "font":{"size": 30,"family": "roboto"}}, 'tickfont': {"size": 20}, 'type' : "log", "range" : [np.log10(log_ymin),np.log10(log_ymax)]}
    #    
    ipl = Plotter()
    ipl.initialise_figure("Current Distribution", x_log_args,y_log_args) 
    #plt.ion()
    # Test the frame is shown properly
    plot_frame(ipl,(0,0.5,1,2,3,20,300,3000),(0,0.5,0.3,0.3,0.3,0.3,0.3,0.3))
    ipl.fig.show()
    time.sleep(0.5)
    #clear_output(wait=True)
    plot_frame(ipl,(0,0.5,1,2,3,20,300,3000),(0,0.5,0.3,0.3,0.1,0.1,0.1,0.1))
    ipl.fig.show()
    time.sleep(0.5)
    #clear_output(wait=True)
    plot_frame(ipl,(0,0.5,1,2,3,20,300,3000),(0,0.2,0.5,0,0.1,0.1,0.1,0.1))
    ipl.fig.show()
    
# # %%
# import plotly.graph_objs as go
# import plotly.io as pio
# import IPython.display as display

# pio.renderers.default = "notebook"

# # Create the plot
# fig = go.Figure()
# fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6], mode='markers'))

# # Show the plot
# plot_html = pio.to_html(fig, full_html=False)
# display_obj = display.display(display.HTML(plot_html), display_id=True)

# # Update the plot
# fig.update_traces(x=[1, 2, 2], y=[4, 2, 6])
# plot_html = pio.to_html(fig, full_html=False)

# # Update the displayed plot
# display.update_display(display.HTML(plot_html), display_id=display_obj.display_id)
# # %%

# import plotly.graph_objs as go
# from plotly.subplots import make_subplots
# from IPython.display import clear_output
# import time

# # Create a new subplot
# fig = make_subplots(rows=1, cols=1)

# # Plot a line chart
# trace = go.Scatter(x=[1, 2, 3], y=[4, 5, 6])
# fig.add_trace(trace)



# # Display the plot
# fig.show()
# time.sleep(0.5)
# # Clear the previous output
# clear_output(wait=True)


# fig.update_traces(x=[1, 1, 1], y=[4, 5, 6])

# fig.show()
# # %%
