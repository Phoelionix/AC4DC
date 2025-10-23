# Assumes no additional continuums tracked other than cascade continuum. 

import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    #thaumatin thing
    'axes.titlesize':8,     # fontsize of the axes title
    'axes.labelsize':8,    # fontsize of the x and y labels   
    'ytick.labelsize':8,
    'xtick.labelsize':8,
    'legend.fontsize':8,    
    'lines.linewidth':1.5,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook
from matplotlib.ticker import (AutoMinorLocator)




####
ELECTRON_DENSITY = True # Whether to use electron density for free distribution plots. Electron energy density if False
###
PLOT_ELEMENT_CHARGE= False #
PLOT_FREE_CONTINUUM = False
PLOT_SPLIT_FREE_CONTINUUMS = False
PLOT_COMBINED_SPLIT_FREE_CONTINUUMS = True
PLOT_FREE_SLICES=False
PLOT_ION_RATIOS=False
PLOT_ION_RATIOS_BARS= False
PLOT_ORBITAL_DENSITIES = False #
PLOT_PHOTO_RATES = False

TIGHT_LAYOUT = False
HORIZONTAL_MODE = False
VERTICAL_MODE = False
TWO_VERTICAL_MODE = False

AUTO_VIEW_PLOT = True

#figures_ext = ".eps"
figures_ext = ".png" #.png


###

#COLUMNWIDTH = 3.34646/2
#COLUMNWIDTH = 3.34646/2
#COLUMNWIDTH = 3.34646*1.35
COLUMNWIDTH = 5.90551/3
#COLUMNWIDTH = 7.24409436834/4
#3.34646
#3.4975
#7.24409436834


#COLUMNWIDTH = 7.24409436834/7 # column of document


COLWIDTH = COLUMNWIDTH*1.1 # column of figure
#ROWHEIGHT = COLWIDTH*16/16 * 3


ROWHEIGHT = COLUMNWIDTH*14/16
#ROWHEIGHT = COLUMNWIDTH*12/16
#ROWHEIGHT = COLUMNWIDTH*14/16
#ROWHEIGHT = COLUMNWIDTH*15/16
#ROWHEIGHT = COLUMNWIDTH*8.5/16

#TODO plotter should have general test to see if a subplot is on edge. And if not don't put axis there (if axes all same lim).

#DPI = 100
DPI = 800
##
END_T = None #-17.8 # None # load data up to this time point (fs). None -> final time in data
##
def main():
    set_highlighted_excepthook()



    # Basic num arguments check
    if  len(sys.argv) < 2:
        print("Usage: python3 generate_plots.py Carbon_1")
        print({"Pass extra arguments to plot for each"})
        exit()
        
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    valid_folder_names= True
    for i, data_folder in enumerate(sys.argv[1:]):
        if not path.isdir(molecular_path+data_folder):
            valid_folder_names = False
            print("\033[91mInput error\033[0m (argument \033[91m"+str(i)+ "\033[0m): folder name not found.")
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    for data_folder in sys.argv[1:]:
        label = data_folder +'_Plt'
        make_some_plots(data_folder,molecular_path,label,dname_Figures,PLOT_ELEMENT_CHARGE,PLOT_ION_RATIOS,PLOT_FREE_CONTINUUM,PLOT_FREE_SLICES,PLOT_ION_RATIOS_BARS,PLOT_ORBITAL_DENSITIES,PLOT_PHOTO_RATES,PLOT_SPLIT_FREE_CONTINUUMS,PLOT_COMBINED_SPLIT_FREE_CONTINUUMS)

def make_some_plots(mol_name,sim_output_parent_dir, label,figure_output_dir, tot_charge=False,bound_ionisation=False,free=False,free_slices=False,bound_ionisation_bar=False,orbital_densities_bar=False,photo_rates = False,split_free=False,combined_split_free=False):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
 # extra homeless options
    load_specific_atoms = None# ["C","Fe_singleShell"]#["C","N","O","S","Gd_fast","Cl","Na"]#["Gd_fast"] #["C"]#["Fe_singleShell","C"] #None #["C","N","O"]
    split_continuums_to_load = "DEFAULT" #"all" #["Gd_fast"] #["C"]#["Fe_singleShell","C"]  # None is not valid. Use: [] 


    ############
    # File/directory names
    #######  
    fname_tot_charge = "tot_charge"
    fname_free = "free"
    fname_HR_style = "HR_style"
    fname_bound_dynamics = "bound_dynamics"

    
    

    
    ##############

    
    if load_specific_atoms is not None:
        label+="_"
        for elem in load_specific_atoms:
            label+=elem
    pl = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY,end_t = END_T,
                 num_subplots=1,split_continuums_to_load="all")
    pl.get_atoms(load_specific_atoms,split_continuums_to_load) # so that plotter knows num continuums
    ### single cascade
    scale_factor = 10
    #energy = 500
    #ymax=energy*6/5
    ymax = 9000
    
    cascade_start_time = -7.5

    # free continuum 
    ax,ax2 = pl.plot_free(log=True,ylog=False,cmin=10**(-8-scale_factor),cmax=10**(-3.609-scale_factor),ylim=[-10,-9], xlim=[cascade_start_time,None],
                    keV=True,every=1,every_e=1,
                    show_cbar=False,
                    show_pulse_profile=True,
                    )
    
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    #ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.locator_params(axis='y',nbins=5)
    ax.locator_params(axis='x',nbins=6)

    # intensity-averaged electrons freed
    show_electrons_freed = False
    if show_electrons_freed:
        pl.num_plotted = 0
        ax3 = ax.twinx()

        c='tab:green'
        #ax3.spines['right'].set_color(c)
        ax3.yaxis.label.set_color(c)
        ax3.tick_params(axis='y',which='both', colors=c)
        #pl.plot_tot_charge(ylim=[None,None],every=1,charge_difference=False,legend_loc="best",intensity_averaged=True,custom_ax =ax3) 
        hline_kwargs = dict(
            color = c,
            linestyle = "dashed",
            alpha=0.5,
        )

        pl.plot_electrons_freed(ylim=[None,230],xlim=[cascade_start_time,None],every=1,custom_ax1_ax2 =(ax3,ax2),normalize_to_time=cascade_start_time,
                                hline_kwargs=hline_kwargs,
                                alpha=0.5,color=c,custom_ylabel = "Num. electrons",
                                show_intensity_averaged=True) 
        ax3.yaxis.tick_right()
        ax3.yaxis.set_label_position("right")

        for a in ax, ax3:
            a.tick_params(which='major',direction='out',pad=1.5,length=2.5)
            a.tick_params(which='minor',direction='out',pad=1.5,length=1.5)
        ax3.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax3.locator_params(axis='y',nbins=8)


    pl.delete_remaining_axes()

    plt.gcf().set_figwidth(COLWIDTH*pl.axs.shape[1])
    plt.gcf().set_figheight(ROWHEIGHT*pl.axs.shape[0])

    
    if TIGHT_LAYOUT:
        plt.tight_layout()
    else:
        #This was a bad idea. You'll likely need to adjust these, sorry.

        #plt.savefig(figure_output_dir + label + figures_ext,dpi=DPI,bbox_inches='tight')
        #pl.fig.subplots_adjust(left=0.135, bottom=0.23, right=0.815, top=0.995) # single continuum
        pl.fig.subplots_adjust(left=0.135, bottom=0.8, right=0.815, top=0.995) # mini intensity plot

    

    plt.savefig(figure_output_dir + label + figures_ext,dpi=DPI,format=figures_ext[1:])
    plt.close()

    if AUTO_VIEW_PLOT:
        # from PIL import Image                                                                                    
        # img = Image.open(figure_output_dir + label + figures_ext)
        # img.show() 
        import os
        os.system("wslview " + figure_output_dir + label + figures_ext)

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option