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


####
ELECTRON_DENSITY = True # Whether to use electron density for free distribution plots. Electron energy density if False
###
PLOT_ELEMENT_CHARGE= True #
PLOT_FREE_CONTINUUM = False
PLOT_SPLIT_FREE_CONTINUUMS = False
PLOT_COMBINED_SPLIT_FREE_CONTINUUMS = False
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
COLUMNWIDTH = 3.484252*0.85
#COLUMNWIDTH = 3.34646*1.35
#COLUMNWIDTH = 7.24409436834/3
#COLUMNWIDTH = 7.24409436834/4
#3.34646
#3.4975
#7.24409436834


#COLUMNWIDTH = 7.24409436834/7 # column of document


COLWIDTH = COLUMNWIDTH # column of figure
#ROWHEIGHT = COLWIDTH*16/16 * 3


#ROWHEIGHT = COLUMNWIDTH*12/16
#ROWHEIGHT = COLUMNWIDTH*14/16
ROWHEIGHT =  max(COLUMNWIDTH*3/4,2.2)*0.8
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
    load_specific_atoms = ["C","N","O"]# ["C","Fe_singleShell"]#["C","N","O","S","Gd_fast","Cl","Na"]#["Gd_fast"] #["C"]#["Fe_singleShell","C"] #None #["C","N","O"]
    split_continuums_to_load = "DEFAULT" #"all" #["Gd_fast"] #["C"]#["Fe_singleShell","C"]  # None is not valid. Use: [] 


    ############
    # File/directory names
    #######  
    fname_tot_charge = "tot_charge"
    fname_free = "free"
    fname_HR_style = "HR_style"
    fname_bound_dynamics = "bound_dynamics"

    
    

    
    ##############
    if split_continuums_to_load == "DEFAULT":
        split_continuums_to_load = "full"
        if combined_split_free or split_free:
            split_continuums_to_load = "all"


    if combined_split_free or split_free:
        assert split_continuums_to_load == "all" or len(split_continuums_to_load)  > 0    
    
    if load_specific_atoms is not None:
        label+="_"
        for elem in load_specific_atoms:
            label+=elem
    pl = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY,end_t = END_T,
                 initialise=False)
    pl.get_atoms(load_specific_atoms,split_continuums_to_load) # so that plotter knows num continuums

    num_atoms = len(pl.atomdict)
    num_subplots = tot_charge + free_slices + bound_ionisation_bar + (bound_ionisation+ orbital_densities_bar+ photo_rates)*num_atoms + free + combined_split_free + split_free*(pl.num_continuums()-1)
    assert(num_subplots > 0) , f"number of sublots is {num_subplots}{', have any plots been specified to plot?' if num_subplots == 0 else ''}"


    assert(VERTICAL_MODE+TWO_VERTICAL_MODE+HORIZONTAL_MODE)<2

    def tmp_func():
        if VERTICAL_MODE:
            plt.cla()
            pl.setup_vertical_axes(num_subplots)
        if TWO_VERTICAL_MODE:
            plt.cla()
            pl.setup_vertical_axes(num_subplots,2)
        if HORIZONTAL_MODE:
            plt.cla()
            pl.setup_horizontal_axes(num_subplots)
    if free:
        pl.initialise(None,num_subplots,"full") # Load full continuum. 
        tmp_func()
        pl.plot_free(log=True,ylog=False,cmin=10**(-8),cmax=10**(-3.609),ylim=[0,8000],keV=True,show_title=False,
                     show_cbar = (split_free==False))
        pl.update_inputs(load_specific_atoms=load_specific_atoms,split_continuums_to_load=split_continuums_to_load)
    else:
        pl.initialise(load_specific_atoms,num_subplots,split_continuums_to_load,)
        tmp_func()

    if num_subplots > 1:
        pl.fig.tight_layout()
        pl.fig.subplots_adjust(left=0.12/pl.axs.shape[0], bottom=None, right=None, top=None, wspace=0.2, hspace=0.4)#hspace=None)

    if tot_charge: 
        #NOTE ensure load_specific_atoms is None or does not exclude `atoms` if `atoms` is passed.
        #pl.plot_tot_charge(ylim=[None,None],every=1,charge_difference=False,legend_loc="best")  
        legend_kwargs = dict(
            borderpad=0.2,handletextpad=0.5,handlelength=1,columnspacing=0.35,
        )
        abdullah_colors =  [
                    "#7d7d7d",
                    "#0050a3",
                    "#eb471e",
                    "#437c34",
                    
                ] 
        pl.plot_tot_charge(colours=abdullah_colors,ylim=[0,3],every=1,charge_difference=False,legend_loc="best",legend_frame=False,profile_height_factor=0.87,legend_kwargs=legend_kwargs,right_aligned=True)  
        #pl.plot_tot_charge(ylim=[None,None],xlim=[-18,18],every=1,charge_difference=False,legend_loc="best")  
        #pl.plot_tot_charge(ylim=[None,None],every=1,charge_difference=True,legend_loc="best",atoms=["C","N","O"])  
        #pl.plot_tot_charge(ylim=[0,6],every=1,charge_difference=False,legend_loc="best",atoms=["C","N","O"])  
        #pl.plot_tot_charge(ylim=[0,6],plot_legend=False,every=1,charge_difference=True,legend_loc="best")  
        #pl.plot_tot_charge(ylim=[0,6],plot_legend=False,every=1,charge_difference=True,legend_loc="best",atoms=["C","N","O","S","Gd_fast"])  
        #pl.plot_tot_charge(ylim=[0,1.1],plot_legend=False,every=1,charge_difference=True,legend_loc="best",atoms=["C","N","O","S","Gd_fast"])  
        #pl.plot_tot_charge(ylim=[0,1.1],plot_legend=False,every=1,charge_difference=True,legend_loc="best",atoms=["C","N","O","S"])  
        #pl.plot_tot_charge(ylim=[0,6],plot_legend=False,every=1,charge_difference=False,legend_loc="best",atoms=["C","N","O","S"])  
        #pl.plot_tot_charge(ylim=[0,2.5],plot_legend=False,every=1,charge_difference=True,scale_intensity=0.939)  #TODO automatically set to charge_difference to True if starting with ions...
        #pl.plot_tot_charge(ylim=[0,3.5],xlim=[-15.5,0.5],plot_legend=False,every=1,charge_difference=True)  #TODO automatically set to charge_difference to True if starting with ions...
 

    if bound_ionisation_bar:
        pl.plot_charges_bar("Cr_LDA",show_pulse_profile=True)
        #plt.gcf().set_figwidth(15)        
    if orbital_densities_bar:
        pl.plot_orbitals_bar(atoms=None,atoms_excluded=None,show_pulse_profile=False,normalise = True,
                             show_cbar_label=True, show_title=False,add_element_to_ylabel=True,
                             show_cbar=True,show_cbar_last=True,
                             show_ylabel = False)
        #pl.plot_orbitals_bar(atoms=None,atoms_excluded=["N","O"],show_pulse_profile=False,normalise = True)
        #pl.plot_orbitals_bar("Gd_fast",show_pulse_profile=True,orbitals=["3p","4p","5p"])
    if photo_rates:
        pl.plot_photoionisation(atoms=None,show_pulse_profile=True)
    if bound_ionisation:
        pl.plot_all_charges(show_pulse_profile=False,ylim=[0,1])

        #Abdallah
        #pl.plot_all_charges(show_pulse_profile=False,xlim=[-40,0],ylim=[0,1])
        #pl.fig.set_figwidth(6)
        #pl.fig.set_figheight(4)
        # #Royle B
        # pl.plot_charges_royle_style("Al", True)
        # pl.fig.subplots_adjust(left=0.2,bottom=0.18,top=0.95)
        # pl.fig.set_figheight(2.5) # total dimensions, for when other plots turned off...
        # pl.fig.set_figwidth(6)
        # # leonov
        # pl.plot_charges_leonov_style("Si",show_pulse_profile=True,xlim=[-20,20],ylim=[0,1.099])
        # pl.fig.set_figwidth(6.662*0.7)  
        # pl.fig.set_figheight(6*0.7)  
    every_e = 1 # every nth energy plotted.
    every_t = 1

    if split_free:
        #pl.plot_free(log=True,cmin=10**(-7.609),cmax=1e-3,ylim=[10,8000])
        photo_then_auger = True
        
        order = list(range(pl.num_continuums()-1))
        if photo_then_auger:
            order = order[::2] + order[1::2]
        for _c in order:
            ylim = [0,8000]
            if pl.get_element_and_e_type(_c)[1]=="Auger":
                ylim = [0,3000]
            print(f"Plotting continuum {_c+1}/{pl.num_continuums()-1} {pl.get_element_and_e_type(_c)}")
            pl.plot_free(log=True,ylog=False,cmin=10**(-8),cmax=10**(-3.609),ylim=[0,8000],
                         keV=True,continuum=_c,every=every_t,every_e=every_e,
                         show_cbar=False,show_time_axis_label=False)
        print("Rendering may take some time...")


        # #Leonov
        # ymax = 9e3
        # pl.plot_free(log=True, cmin=10**(-6.609),cmax = 10**(-2), every=5,mask_below_min=True,cmap='turbo',ymax=ymax,leonov_style=True)
        # pl.fig.set_figwidth(6.662*0.7*1.16548042705)  
        # pl.fig.set_figheight(6*0.7)  
        
    if combined_split_free:
        ### full continuum (if want to confirm same as combining all)
        pl.update_inputs(load_specific_atoms=load_specific_atoms,split_continuums_to_load='full')
        pl.plot_free(log=True,ylog=False,cmin=10**(-8),cmax=10**(-3.609),ylim=[0,8000],keV=True,show_title=True,
                     show_cbar = True)
        ### Combined
        # pl.plot_free(log=True,ylog=False,cmin=10**(-8),cmax=10**(-3.609),ylim=[0,8000],
        #              keV=True,every=every_t,every_e=every_e,
        #              show_cbar=False)
    

    if free_slices:
        normed=False

        pl.initialise_step_slices_ax()
        from plotter_core import fit_maxwell, maxwell

        cmap = plt.get_cmap("tab10")

        plt.rcParams["font.size"] = 10 # 8

        ### defaults
        #Cutoff energies for fitting MB curves. TODO get from AC4DC
        thermal_cutoff_energies = [2000]*10     
        xmin,xmax = 1,1e4
        
        ###

        # slices = [-15,0,15]
        #TODO get slices from AC4DC or user input           
        # # Hau-Riege (Sanders results)
        # thermal_cutoff_energies = [200, 500, 500, 1000]
        # slices = [-7.5,-5,-2.5,0]         
        # # -7.5 fs Hau-Riege
        # thermal_cutoff_energies = [200]
        # slices = [-7.5] 

        #abdallah
        #slices = [-39,-38,-36,-34,-32,-30,-0.01]         
        # Royle Sect. B
        #thermal_cutoff_energies = [500,1000,1250,2000]
        #slices = [-90, -50, 0, 50]  
        # # Royle Sect. C
        thermal_cutoff_energies = [500,500,600,2000]
        slices = [-15, 0, 15, 30]
        
        # thermal_cutoff_energies = [2000,2000,2000,2000]
        # slices = [30,65,100]

        colrs = [cmap(i) for i in range(len(slices))]

        plot_legend = True
        plot_fits = False # Whether to fit MB curves to distribution below thermal cutoff energies.
        plot_those_darn_knots = False
        ####### 
        # Here we can plot MB curves e.g. for fit comparison
        plot_custom_fits = False
        ####
        # # example
        # custom_T = [30.8,75.5,131.8,207.5]
        # custom_n = [0.06*3/2,0.06*3/2,0.06*3/2,0.06*3/2]
        # #Sanders/H-R for -7.5 fs  
        # custom_T = [44.1,31]  
        # custom_n = [0.06*3/2,0.06*3/2] # - (not anything meaningful, just density of MB approximated as 50% of total). 
        # custom_colrs = ['r','b']
        #H-R
        # custom_T = [31,70,125,195]   # [44.1,84.9,135.6,205.8]
        # custom_n = [0.06*3/2]*len(custom_T) # - (not anything meaningful, just density of MB approximated as 50% of total). 
        # custom_colrs = colrs   
        # xmin,xmax = 1,1e4 
        # plot_legend = True 
        #Royle B
        # xmin, xmax = 0,2000
        # colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]  
        # custom_colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]  
        # plot_fits = True
        # plot_custom_fits = False
        # custom_T = [8,33,105,120]  
        # custom_n = [0.12,0.05,0.12,0.12] 
        # #Royle C
        xmin,xmax = 0,1e4
        custom_T = [100]  
        custom_n = [0.155]     
        colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]          
        custom_colrs = ['black']
        plot_fits = False
        plot_custom_fits = False
        # Fit our last point in time 
        if plot_custom_fits == False:
            T = pl.plot_fit(slices[-1], thermal_cutoff_energies[-1], normed=normed, color=cmap(3), lw=1.5,alpha=0.7)
        #######

        #v_anchors = [0.16,0.12,0.08,0.04]
        #v_anchors = [0.2,0.15,0.1,0.05]
        v_anchors = [0.19,0.12,0.05]

        lw = 1.5
        #pl.ax_steps.set_ylim([0.4e-4, 0.4])
        pl.ax_steps.set_ylim([1e-4*1e3, 1*1e3])
        #pl.ax_steps.set_ylim([1e-4, 1])
        #TODO get from AC4DC
        pl.ax_steps.set_xlim([xmin,xmax]) #Hau-Riege        
        lines = []
        for (t, e, col ) in zip(slices, thermal_cutoff_energies, colrs):
            lines.extend(pl.plot_step(t, normed=normed, color = col, lw=lw))
            if plot_fits:
                T = pl.plot_fit(t, e, normed=normed, color=col, lw=lw,alpha=0.7)



        if plot_custom_fits:
            for (T, n, col ) in zip(custom_T, custom_n, custom_colrs):
                pl.ax_steps.plot([0],[0],alpha=0,label=None)
                pl.plot_maxwell(T,n,color = col, lw=lw,alpha=0.7)
        if plot_those_darn_knots:
            ymin,_ = pl.ax_steps.get_ylim()
            pl.ax_steps.set_ylim([ymin*0.67,None])
            pl.plot_the_knots(slices,v_anchors,colrs,padding=0.14)
            #pl.fig_steps.set_figheight(4.8)
            # pl.fig_steps.set_figwidth(7)  
            # for l in lines:
            #     l.set_linewidth(3)          
        
        if ELECTRON_DENSITY:
            #pl.ax_steps.set_ylim([2e-7*1e3, 1e-2*1e3]) #royle sect. B
            pl.ax_steps.set_ylim([1e-8*0.8e3, 1e-3*0.8e3]) #royle sect. C


        #pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.95)
        pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
        pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False


        if plot_legend:
            handles, labels = pl.ax_steps.get_legend_handles_labels()
            ncols = 2
            # e.g. for 8 labels (4 time steps), order = [0,2,4,6,1,3,5,7]  , if ncols = 2.
            order = []
            #2 cols
            order = list(range(0,len(labels) - 1,ncols)) + list(range(1,len(labels),ncols)) 
            # 2 rows 4 cols
            #order = list(range(0,len(labels) - 1,ncols)) + list(range(1,len(labels),ncols)) +  list(range(2,len(labels) - 1,ncols)) + list(range(ncols-1,len(labels),ncols)) 
            if len(labels)%2 != 0:
                order.append(len(order))  # shouldnt happen though.
            if plot_fits == False and plot_custom_fits == False:
                ncols = 1
                order = list(range(0,len(labels)))
            leg = pl.ax_steps.legend([handles[idx] for idx in order],[labels[idx] for idx in order],ncol=ncols,
                                     loc='upper right',bbox_to_anchor=(0.85, 1),borderpad=0,labelspacing =0.1)
            leg.get_frame().set_linewidth(0)
            

        name = label.replace('_',' ')
        pl.ax_steps.set_title(name + " - Free-electron distribution")

        #plt.savefig(figure_output_dir + label + fname_HR_style + figures_ext)
    #Abdallah
    # pl.ax_steps.set_xscale("linear")
    # pl.ax_steps.set_xlim([0,1800])
    # pl.ax_steps.set_ylim([0.5e-4,0.5])
    pl.delete_remaining_axes()

    plt.gcf().set_figwidth(COLWIDTH*pl.axs.shape[1])
    plt.gcf().set_figheight(ROWHEIGHT*pl.axs.shape[0])

    
    if TIGHT_LAYOUT:
        plt.tight_layout()
    else:
        #This was a bad idea. You'll likely need to adjust these, sorry.

        #plt.savefig(figure_output_dir + label + figures_ext,dpi=DPI,bbox_inches='tight')
        #pl.fig.subplots_adjust(left=0.185, bottom=0.22, right=0.975, top=0.945, wspace=0.2, hspace=0.4) # Custom 
        pl.fig.subplots_adjust(left=0.01,bottom=0.2,top=0.9,right=0.87)
        if VERTICAL_MODE:
            # keeps canvas height (almost) the same for 1 or 3 subplots 
            plt.gcf().set_figheight(ROWHEIGHT*pl.axs.shape[0]*(0.9+0.205/num_subplots)) 
            pl.fig.subplots_adjust(left=0.24, bottom=0.42/num_subplots, right=0.95, top=1-0.005/num_subplots, wspace=0.2, hspace=0.4) # Custom for column of orbital density plots 
        if HORIZONTAL_MODE:
            # keeps canvas height (almost) the same for 1 or 3 subplots 
            #plt.gcf().set_figwidth(ROWHEIGHT*1.1) 

            #colwidth 11/16 row width, row width 3.5/2
            # plt.gcf().set_figwidth(COLWIDTH*pl.axs.shape[1]*(0.9+0.205/num_subplots))
            # pl.fig.subplots_adjust(left=0.25*3**1.2/3/num_subplots**1.2, bottom=0.42, right=1-0.09/num_subplots, top=0.998, wspace=0.43, hspace=0.4) # Custom for row of orbital density plots 
            
            # Don't bother, just plot dummy plots to make it fit...
            #plt.gcf().set_figwidth(COLWIDTH*pl.axs.shape[1]*(0.9+0.3/num_subplots))
            #pl.fig.subplots_adjust(left=0.21*3**1.2/3/num_subplots**1.2, bottom=0.41, right=1-0.005/num_subplots, top=0.996, wspace=0.38, hspace=0.4) # Custom for row of orbital density plots 
            #pl.fig.subplots_adjust(left=0.025, bottom=0.415, right=0.7, top=0.996, wspace=0.3, hspace=0.4) # 3 orb density plots row
            pl.fig.subplots_adjust(left=0.08, bottom=0.415, right=0.83, top=0.93, wspace=0.3, hspace=0.4) # 3 orb density plots row
        
        if TWO_VERTICAL_MODE:
            pl.fig.subplots_adjust(left=0.1, bottom=0.03, right=0.95, top=0.97, wspace=0.16, hspace=0.6) # Custom for column of orbital density plots 
    
    

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