#%% 

#sys_argv = ["","Carbon_165_50"]
sys_argv = ["","Sulfur_Sanders_Settings"]

import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
#%matplotlib widget
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from plotter_core import Plotter, fit_maxwell, maxwell
import numpy as np
import pickle
import sys, traceback

# Exception colouring
# https://stackoverflow.com/questions/14775916/coloring-exceptions-from-python-on-a-terminal
def set_highlighted_excepthook():  
    from pygments import highlight
    from pygments.lexers import get_lexer_by_name
    from pygments.formatters import TerminalFormatter

    lexer = get_lexer_by_name("pytb" if sys.version_info.major < 3 else "py3tb")
    formatter = TerminalFormatter()

    def myexcepthook(type, value, tb):
        tbtext = ''.join(traceback.format_exception(type, value, tb))
        sys.stderr.write(highlight(tbtext, lexer, formatter))

    sys.excepthook = myexcepthook
set_highlighted_excepthook()

# Basic num arguments check
if  len(sys_argv) < 2:
    #print("Usage: python3 _save_interactive.py Carbon_1 <number of steps> <name_suffix (optional)>")
    print("Usage: python3 _save_interactive.py Carbon_1 <name_suffix (optional)>")
    exit()


label = sys_argv[1]
if len(sys_argv) > 2:
     label += '_' + sys_argv[2]
name = sys_argv[1].replace('_',' ')
pl = Plotter(sys_argv[1],"y")

#ParamsthatshouldreallygoinafilesomewhereO_0
pl.yscale = "linear"
outdir = "../../AC4DC_Interactives/"
pl.lin_ymax = 0.035
pl.lin_ymin = -0.035
pl.log_ymax = 1
pl.log_ymin = 1e-10
pl.xmax = 1e4
pl.xmin = 1
normalise = True

# Initialisation
pl.max_t = pl.timeData[-1]
pl.min_t = pl.timeData[0]
pl.fig_steps.set_size_inches(6,5)
pl.ax_steps.set_xlim([pl.xmin,pl.xmax]) 
pl.line = pl.plot_step(pl.min_t, normed=normalise, lw=0.7, color=(0,0,0,1))
#nsteps = int(sys_argv[2])
#times = np.linspace(pl.min_t,pl.max_t,nsteps)

#-----Widgets-----#
# Time Slider
pl.ax_steps.set_yscale('linear')
scale = 5
pl.axtime = pl.fig_steps.add_axes([0.25, 0.1, 0.65*scale, 0.03*scale])
time_slider = Slider(
    ax=pl.axtime,
    label='Time [fs]',
    valmin=pl.min_t,
    valmax=pl.max_t,
    valinit=pl.min_t,
)

# Y-Scale Button
pl.scale_button_ax = pl.fig_steps.add_axes([0.8, 0.025, scale*0.1, scale*0.04])
scale_button = RadioButtons(pl.scale_button_ax, ['Log','Lin','SymLog'],0, activecolor='0.975') ## ミ=͟͟͞͞(✿ʘ ᴗʘ)っ

def widget_update(val):  
    t = time_slider.val
    ### Cute colours  ★ﾟ~(◠ᴗ ◕✿)ﾟ*｡:ﾟ+ 
    rgb_intensity = [1/2,2/3,1]  # max = 1
    rgb_width = [1,1,0.5]
    rgb = [0,0,1]
    rgb_direction = [1,1,-1]
    # Linear interpolation
    rgb_bndry = [1/x for x in rgb_width]
    rgb[0] =  (1/rgb_width[0])*((t-pl.min_t)/(pl.max_t-pl.min_t)-rgb_bndry[0])
    rgb[1] =  (1/rgb_width[1])*((t-pl.min_t)/(pl.max_t-pl.min_t)-rgb_bndry[1])
    rgb[2] =  (1/rgb_width[1])*((t-pl.min_t)/(pl.max_t-pl.min_t)-rgb_bndry[1])
    for i in range(len(rgb)): 
        if rgb_direction[i] < 0:
            rgb[i] = 1 - rgb[i]
        rgb[i] = min(1, max(0, rgb[i])*rgb_intensity[i]) 
    ### Update line
    pl.line[-1].remove()
    col = tuple(rgb)
    pl.ax_steps.set_ylim([pl.log_ymin, pl.log_ymax])
    pl.line = pl.plot_step(t, normed=normalise, lw=0.7, color=col)     # pl handles the continuous values of the time for us. Also sets the x-y scales to log.
    if scale_button.value_selected == "Log":
        pass
    if scale_button.value_selected == 'Lin':
        pl.ax_steps.set_yscale('linear')     
        pl.ax_steps.set_ylim([pl.lin_ymin, pl.lin_ymax])
    if scale_button.value_selected == 'SymLog':
        pl.ax_steps.set_yscale('symlog',linthresh = 1e-10)
        pl.ax_steps.set_ylim([-pl.log_ymax, pl.log_ymax])        
    pl.fig_steps.canvas.draw_idle()  
pl.widget_update = widget_update

pl.widget_update(0) 
scale_button.on_clicked(pl.widget_update)
time_slider.on_changed(pl.widget_update)
time_slider.set_val(pl.min_t+(pl.max_t-pl.min_t)/10)
scale_button.set_active(0)


#-----Plot-----#

pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.95)
pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False

#pl.plot_maxwell(44.1,0.06*3/2)  # ~Sanders -7.5 fs - density of MB assumed to be 50% of total.  

#handles, labels = pl.ax_steps.get_legend_handles_labels()

#pl.ax_steps.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper left',ncol=2)
#plt.gcf()

plt.rcParams["font.size"] = 16 # 8
plt.title(name + " - Free-electron distribution")
#plt.tight_layout()

#plt.savefig(outdir + "test")
extension = ".fig.pickle"

#pickle.dump(pl, open(outdir + label + extension, 'wb'))

#plotdata = {'pl':pl,'scale_button':scale_button,'time_slider':time_slider}
#pickle.dump(plotdata, open(outdir + label + extension, 'wb'))

print("Done!")

#TODO: Change size to match my screen by default, add --y option
# %%
