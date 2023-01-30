#%% 


### Params

#sys_argv = ["","Carbon_165_50"]
sys_argv = ["","Sulfur_Sanders_Settings"]
xmax = 1e4


# import matplotlib
# matplotlib.use("pgf")
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })
#matplotlib widget
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from plotter_core import Plotter #, fit_maxwell, maxwell
import pickle
import sys, traceback

# Exception colouring
# https://stackoverflow.com/questions/14775916/coloring-exceptions-from-python-on-a-terminal
def set_highlighted_excepthook_c():  
    from pygments import highlight
    from pygments.lexers import get_lexer_by_name
    from pygments.formatters import TerminalFormatter

    lexer = get_lexer_by_name("pytb" if sys.version_info.major < 3 else "py3tb")
    formatter = TerminalFormatter()

    def myexcepthook(type, value, tb):
        tbtext = ''.join(traceback.format_exception(type, value, tb))
        sys.stderr.write(highlight(tbtext, lexer, formatter))

    sys.excepthook = myexcepthook
set_highlighted_excepthook_c()

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

outdir = "../../AC4DC_Interactives/"

normalise = True

# Initialisation
pl.max_t = pl.timeData[-1]
pl.min_t = pl.timeData[0]
pl.fig_steps.set_size_inches(6,5)
pl.ax_steps.set_xlim([1,xmax]) 
pl.line = pl.plot_step(pl.min_t, normed=normalise, lw=0.7, color=(0,0,0,1))
#nsteps = int(sys_argv[2])
#times = np.linspace(pl.min_t,pl.max_t,nsteps)

#-----Widgets-----#
# Time Slider
pl.ax_steps.set_yscale('linear')
scale = 5
pl.axtime = pl.fig_steps.add_axes([0.25, 0.1, 0.65*scale, 0.03*scale])
# time_slider = Slider(
#     ax=pl.axtime,
#     label='Time [fs]',
#     valmin=pl.min_t,
#     valmax=pl.max_t,
#     valinit=pl.min_t,
# )

# # Y-Scale Button
pl.scale_button_ax = pl.fig_steps.add_axes([0.8, 0.025, scale*0.1, scale*0.04])
# scale_button = RadioButtons(pl.scale_button_ax, ['Log','Lin','SymLog'],0, activecolor='0.975') ## ミ=͟͟͞͞(✿ʘ ᴗʘ)っ

class WidgetHandler_c:  
    def __init__(self,_pl,_time_slider = None,_scale_button = None):
        #ConstantsthatshouldreallygoinafilesomewhereandmaybecapitalisedO_0
        self.lin_ymax = 0.035
        self.lin_ymin = -0.035
        self.log_ymax = 1
        self.log_ymin = 1e-10      
        #-----Widgets-----#
        # Time Slider        
        if _time_slider == None:
            _pl.ax_steps.set_yscale('linear')
            _time_slider = Slider(ax=_pl.axtime, 
            label='Time [fs]', 
            valmin=_pl.min_t, 
            valmax=_pl.max_t, 
            valinit=_pl.min_t,
            )
        
            #time_slider.set_val(pl.min_t+(pl.max_t-pl.min_t)/10)
        if _scale_button == None:
            _scale_button = RadioButtons(_pl.scale_button_ax, 
            ['Log','Lin','Logsym'],
            0,
            activecolor='tab:orange' ## ミ=͟͟͞͞(✿ʘ ᴗʘ)っ
            ) 

        # Connect to widgets
        _time_slider.on_changed(self.update)
        _scale_button.on_clicked(self.update)

        # copy references to sliders and plot objects
        self.time_slider = _time_slider
        self.scale_button = _scale_button
        self._pl = _pl 
    def update(self, dummy_val):
        pl = self._pl
        time_slider = self.time_slider
        scale_button = self.scale_button
        lin_ymax, lin_ymin, log_ymax, log_ymin, = self.lin_ymax, self.lin_ymin, self.log_ymax, self.log_ymin,

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
        pl.ax_steps.set_ylim([log_ymin, log_ymax]) # Needed as pl.plot_step sets scales to log by default.
        pl.line = pl.plot_step(t, normed=normalise, lw=0.7, color=col)     # pl handles the continuous values of the time for us.
        if scale_button.value_selected == "Log":
            pass
        if scale_button.value_selected == 'Lin':
            pl.ax_steps.set_yscale('linear')     
            pl.ax_steps.set_ylim([lin_ymin, lin_ymax])
        if scale_button.value_selected == 'SymLog':
            pl.ax_steps.set_yscale('symlog',linthresh = 1e-10)
            pl.ax_steps.set_ylim([-log_ymax, log_ymax])        
        pl.fig_steps.canvas.draw_idle()  

# widg_test = WidgetHandler(pl,time_slider,scale_button)

# widg_test.update(0) 
# scale_button.on_clicked(widg_test.update)
# time_slider.on_changed(widg_test.update)
# time_slider.set_val(pl.min_t+(pl.max_t-pl.min_t)/10)
# scale_button.set_active(0)


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

pickle.dump(pl, open(outdir + label + extension, 'wb'))
#plotdata = {'pl':pl,'pl.scale_button':pl.scale_button,'time_slider':time_slider}
#pickle.dump(plotdata, open(outdir + label + extension, 'wb'))

print("Done!")
plt.close() # IMPORTANT. If importing methods the whole code is run.

#TODO: Change size to match my screen by default, add --y option
# %%

# %%
