import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
from plotter_core import fit_maxwell, maxwell
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook

ELECTRON_DENSITY = False # energy density if False
NORMED = False

def main():
    set_highlighted_excepthook()

    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    valid_folder_names= True
    for i, data_folder in enumerate(sys.argv[1:]):
        if not path.isdir(molecular_path+data_folder):
            valid_folder_names = False
            print("\033[91mInput error\033[0m (argument \033[91m"+str(i)+ "\033[0m): folder name not found.")
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    data_folders = sys.argv[1:]
    label = data_folders[0] +'..._Plt-slices'
    
    make_the_plot(data_folders,molecular_path,label,dname_Figures)

def make_the_plot(mol_names,sim_output_parent_dir, label,figure_output_dir):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
    ############
    # File/directory names
    #######  
    figures_ext = ".pdf" #.png
    fname_charge_conservation = "charge_conservation"
    fname_free = "free"
    fname_HR_style = "HR_style"
    fname_bound_dynamics = "bound_dynamics"

    cmap = plt.get_cmap("tab10")
    num_slices = 1
    all_colrs = [cmap(i) for i in range(num_slices*len(mol_names))] # (needlessly generalised)

    fig_steps, (ax_steps1,ax_steps2) = plt.subplots(1, 2, sharey=False, facecolor='w')
    for m, mol_name in enumerate(mol_names):           
        colrs = all_colrs[m*num_slices:(m+1)*num_slices]
        pl1 = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY)
        plt.subplots_adjust(hspace=0,wspace=0) 

        pl2 = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY)      
        plt.subplots_adjust(hspace=0,wspace=0) 

        
        pl1.fig_steps, (pl1.ax_steps,pl2.ax_steps) = fig_steps, (ax_steps1,ax_steps2)
        
        pl2.fig_steps = pl1.fig_steps

        
        

        # hide the spines between ax and ax2
        pl1.ax_steps.spines['right'].set_visible(False)
        pl2.ax_steps.spines['left'].set_visible(False)
        pl1.ax_steps.yaxis.tick_left()
        #pl1.ax_steps.tick_params(labelright='off')
        pl2.ax_steps.yaxis.tick_right()    
        #pl2.ax_steps.yaxis.set_visible(False)
        pl2.ax_steps.tick_params(labelright='off')
        # # y version:
        # # pl1.ax_steps.spines['bottom'].set_visible(False)
        # # pl2.ax_steps.spines['top'].set_visible(False)
        # # pl1.ax_steps.xaxis.tick_top()
        # # #pl2.ax_steps.tick_params(labelbottom='off')
        # # pl2.ax_steps.xaxis.tick_bottom()

        plt.rcParams["font.size"] = 10 # 8

        ### defaults
        #Cutoff energies for fitting MB curves. TODO get from AC4DC
        thermal_cutoff_energies = [2000]*10     
        
        # Hau-Riege (Sanders results)
        #thermal_cutoff_energies = [200, 500, 500, 1000]
        #slices = [-7.5,-5,-2.5,-0.01]   
        # xmin1,xmax1 = 1,1e3
        # xmin2,xmax2 = 3e3,1e4    
        thermal_cutoff_energies = [1000]*num_slices
        ## (not normed)
        # slices = [-7.5]
        # xmin1,xmax1 = 1,400
        # xmin2,xmax2 = 4e3,7e3    
        # ylim1 = 0.03; ylim2 = 0.3
        slices = [-9.75]
        xmin1,xmax1 = 1,500
        xmin2,xmax2 = 5.5e3,6e3   
        ylim1 = 0.00175; ylim2 = 0.35
        # slices = [-9]
        # xmin1,xmax1 = 1,400
        # xmin2,xmax2 = 5.2e3,6.2e3   
        # ylim1 = 0.025; ylim2 = 0.58

        # slices = [-7.5]
        # xmin1,xmax1 = 10,1000
        # xmin2,xmax2 = 5.2e3,6.2e3   
        # ylim1 = 0.056; ylim2 = 0.51
        ## normed 
        # slices = [-9.5]
        # xmin1,xmax1 = 1,400
        # xmin2,xmax2 = 5.2e3,6.2e3   
        # ylim1 = 0.024; ylim2 = 1.5      
        # slices = [-7.5]
        # xmin1,xmax1 = 1,400
        # xmin2,xmax2 = 5.2e3,6.2e3   
        # ylim1 = 0.03; ylim2 = 0.3

        pl1.ax_steps.set_xlim([xmin1,xmax1])    
        pl2.ax_steps.set_xlim([xmin2,xmax2])    
        pl1.ax_steps.set_ylim([0, ylim1])
        pl2.ax_steps.set_ylim([0, ylim2])
        
        # Break lines, but looks kinda bad.
        # d = .012  # length of break lines
        # kwargs = dict(transform=pl1.ax_steps.transAxes, color='k', clip_on=False,lw=0.5)
        # pl1.ax_steps.plot((1-d, 1+d), (-d, +d), **kwargs)
        # pl1.ax_steps.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

        # kwargs.update(transform=pl2.ax_steps.transAxes)  # switch to the bottom axes
        # pl2.ax_steps.plot((-d, +d), (1-d, 1+d), **kwargs)
        # pl2.ax_steps.plot((-d, +d), (-d, +d), **kwargs)


        

        plot_legend = False       
        plot_fits = False # Whether to fit MB curves to distribution below thermal cutoff energies.
        ####### 
        # Here we can plot MB curves e.g. for fit comparison
        ####
        plot_custom_fits = False
        
        for i, pl in enumerate((pl1,pl2)):
            lw = 1.5
            for (t, e, col ) in zip(slices, thermal_cutoff_energies, colrs):
                    lines = pl.plot_step(t, normed=NORMED, color = col, lw=lw)
                    if plot_fits:
                        T = pl.plot_fit(t, e, normed=NORMED, color=col, lw=lw,alpha=0.7)
            if plot_custom_fits:
                for (T, n, col ) in zip(custom_T, custom_n, custom_colrs):
                    pl.ax_steps.plot([0],[0],alpha=0,label=" ")
                    pl.plot_maxwell(T,n,color = col, lw=lw,alpha=0.7)


            

            #pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.95)
            pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
            pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False

            if plot_legend:
                handles, labels = pl.ax_steps.get_legend_handles_labels()
                # e.g. for 8 labels (4 time steps), order = [0,2,4,6,1,3,5,7]  
                order = list(range(0,len(labels) - 1,2)) + list(range(1,len(labels),2))
                if len(labels)%2 != 0:
                    order.append(len(order))  # shouldnt happen though.

            
                pl.ax_steps.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper center',ncol=2)

            # name = label.replace('_',' ')
            # pl.ax_steps.set_title(name + " - Free-electron distribution")
            #plt.savefig(figure_output_dir + label + fname_HR_style + figures_ext)
            # pl.ax_steps.set_xscale(xscale)
            #pl.ax_steps.set_yscale('log')            
            pl.ax_steps.set_yscale('linear')
            pl1.ax_steps.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
            pl2.ax_steps.ticklabel_format(axis='y', style='sci',scilimits=(0,0))          
    
    pl.fig_steps.set_figheight(5)
    pl.fig_steps.set_figwidth(16) 
    plt.tight_layout()
    plt.savefig(figure_output_dir + label + figures_ext)
    plt.close()
    

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option

    # pl1 = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY)
    # pl2 = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY)
    
    # pl1.fig_steps, (pl1.ax_steps,pl2.ax_steps) = plt.subplots(1, 2, sharex=True, facecolor='w')
    # pl2.fig_steps = pl1.fig_steps


    # cmap = plt.get_cmap("tab10")

    # plt.rcParams["font.size"] = 10 # 8

    # ### defaults
    # #Cutoff energies for fitting MB curves. TODO get from AC4DC
    # thermal_cutoff_energies = [2000]*10     
    # xmin,xmax = 1,1e4
    
    # # Hau-Riege (Sanders results)
    # thermal_cutoff_energies = [200, 500, 500, 1000]
    # slices = [-7.5,-5,-2.5,-0.01]   

    # colrs = [cmap(i) for i in range(len(slices))]

    # plot_legend = True        
    # plot_fits = False # Whether to fit MB curves to distribution below thermal cutoff energies.
    # ####### 
    # # Here we can plot MB curves e.g. for fit comparison
    # ####
    # plot_custom_fits = False
    
    # for pl in (pl1,pl2):
        