# Plotter object determines electron density f(e), then we multiply by sqrt(e)*sigma_EII(e) 

import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': 14,
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
from labellines import *
from copy import deepcopy 

NORMED = False
FULL_FIG = True
PLOT_RATE_SEPARATELY = False


def main():
    set_highlighted_excepthook()
    assert len(sys.argv[1:]) == 2, "Usage: python EII_rate_plot.py <sim_handle_1> <sim_handle_2>"
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

def make_the_plot(mol_names,sim_output_parent_dir, label,figure_output_dir,plot_rate_separately = False):
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

    
    if FULL_FIG:
        fig_steps, (ax_steps1) = plt.subplots(1, 1, sharey=False, facecolor='w')
    else:
        fig_steps, (ax_steps1,ax_steps2) = plt.subplots(1, 2, sharey=False, facecolor='w')

    dashes = ["dashed","solid",]#"dotted"]
    assert len(mol_names) <= 2

    # CARBON EI CROSS-SECTION
    # Plot the electron-impact ionisation cross-section from Suno&Kato https://doi.org/10.1016/j.adt.2006.01.001
    # fit parameters for each transition ("0,1" means from charge of 0 to charge of 1)
    # alternative but annoying and possibly costly method would be to incorporate rate tracking in the ODE solver.
    transitions =  dict()
    # charge                I      A1        A2        A3        A4        A5      rms 
    transitions["0,1"] =  "10.6 1.829E+0 -1.975E+0  1.149E+0 -3.583E+0  2.451E+0   0.61".split()
    transitions["1,2"] =  "24.4 8.390E-1 -7.950E-1  3.263E+0 -5.382E+0  3.476E+0   0.24".split()  
    for key,val in transitions.items():
        transitions[key] = [float(s) for s in val] 
    # sigma = {}
    # for key, V in transitions.items():
    #     I = V[0]; A = V[1:6]; rms = V[6]
    #     sigma[key] = lambda E: 10e-13/(I*E)*(A[0]*np.log(E/I)+np.sum( [A[i]*(1-I/E)**(i) for i in range(1,len(A))] ))
    def sigma(E,transition):
        V = transitions[transition]
        I = V[0]; A = V[1:6]; rms = V[6]
        if E <= I: 
            return 0
        return 1e-13/(I*E)*(A[0]*np.log(E/I)+np.sum( [A[i]*(1-I/E)**(i) for i in range(1,len(A))] ))
    def eii_prefactor(E,transition):
        ang_per_cm = 1e8
        fs_per_eV = 0.6582
        return sigma(E,transition)*np.sqrt(E)*ang_per_cm**2/fs_per_eV # Convert to cross-section in angstrom and add conversion of rate to fs.  
    
    slices = [-9.5,-8.5,-7.5]
    xmin,xmax = 10,1e4

    electron_density = (PLOT_RATE_SEPARATELY == False)
    label_x_anchors_override = [22,70,90,22,50,72]
    for m, mol_name in enumerate(mol_names):           
        pl1 = Plotter(mol_name,sim_output_parent_dir,use_electron_density = electron_density)
        pl2 = Plotter(mol_name,sim_output_parent_dir,use_electron_density = electron_density)      
       

        pl1.fig_steps, pl1.ax_steps, = fig_steps, ax_steps1
        plotters = (pl1,)   
        if not FULL_FIG:
            pl2.ax_steps = ax_steps2
            pl2.fig_steps = pl1.fig_steps
            plotters += (pl2,)

        ### defaults
        #Cutoff energies for fitting MB curves. TODO get from AC4DC
        thermal_cutoff_energies = [2000]*10     

        if PLOT_RATE_SEPARATELY and m == 0:
            y = []
            E = np.logspace(1,4,2000)
            # can't simply vectorise for whatever reason.
            for e in E:
                y.append(sigma(e,"0,1")/np.sqrt(e))   # f(e)*sqrt(e) = cross-section, so need to divide by sqrt(e) to see effecto f f(e)*e
            y = np.array(y)
            ax_eii = pl1.ax_steps.twinx()
            ax_eii.tick_params(axis='y', colors='blue')
            ax_eii.set_ylabel("$\sigma^{EI}(\epsilon) \epsilon^{-1/2}$",color="blue")
            ax_eii.plot(E[y>=0],y[y>=0],color="black",linestyle="dotted")
            ax_eii.set_ylim(0)


        cmap = plt.get_cmap("tab10")
        colrs = [cmap(i) for i in range(len(slices))] # (needlessly generalised)        
        linestyle = dashes[m]

        thermal_cutoff_energies = [1000]*len(slices)


        plot_legend = True     
        plot_fits = False # Whether to fit MB curves to distribution below thermal cutoff energies.
        ####### 
        # Here we can plot MB curves e.g. for fit comparison
        ####
        plot_custom_fits = False

        prefactor_kwargs = {}
        if not PLOT_RATE_SEPARATELY:
            prefactor_kwargs = dict(prefactor_function = eii_prefactor, prefactor_args = dict(transition="0,1"))
        for i, pl in enumerate(plotters):
            lw = 2
            lines = []
            for (t, e, col ) in zip(slices, thermal_cutoff_energies, colrs):
                    lines.extend(pl.plot_step(t, normed=NORMED, color = col, lw=lw,linestyle=linestyle,**prefactor_kwargs))
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
                
                pl.ax_steps.legend(handles[len(slices):],labels[len(slices):], loc='upper center',ncol=1)
            

            
            energies_vect = [line.get_xdata() for line in lines]
            y_vect = [line.get_ydata() for line in lines]

            integrated_rates = []
            label_x_anchors = []
            for i in range(len(y_vect)):
                energies = energies_vect[i]; y = y_vect[i]
                de = np.append(energies, energies[-1]*2 - energies[-2]) 
                de = de[1:] - de[:-1]
                integrated_rates.append("{:.2f}".format(np.dot(y, de)))
                label_x_anchors.append(energies[np.array(y).argmax()]*3.3)
            
            #ovverriding..
            label_x_anchors = label_x_anchors_override[len(slices)*m:len(slices)*(m+1)]
            labelLines(lines,xvals=label_x_anchors,forced_labels=integrated_rates)  #  ATTENTION: I was lazy and hacked this module with forced_labels  

            # name = label.replace('_',' ')
            # pl.ax_steps.set_title(name + " - Free-electron distribution")
            #plt.savefig(figure_output_dir + label + fname_HR_style + figures_ext)
            # pl.ax_steps.set_xscale(xscale)
            #pl.ax_steps.set_yscale('log')            
            pl.ax_steps.set_yscale('linear')
            pl.ax_steps.ticklabel_format(axis='y', style='sci',scilimits=(0,0))

            pl.ax_steps.set_ylim([0,None])
            pl.ax_steps.set_xlim([xmin,xmax])
            #pl.ax_steps.set_ylabel("$\epsilon^{1/2} \sigma(\epsilon) f(\epsilon)$",rotation='horizontal')
            #pl.ax_steps.set_ylabel("$\\frac{d}{d\epsilon}\Gamma^{EI}(\epsilon)$",rotation='horizontal')
            #pl.ax_steps.set_ylabel("Impact ionisation rate $\\left(\\frac{d\\,\\Gamma^{ei}}{d\\epsilon}\\right)$")
            #pl.ax_steps.set_ylabel("Impact ionisation rate ($\\textrm{fs}^{-1}$)")
            pl.ax_steps.set_ylabel("Impact ionisation rate (fs$^{-1}$)")
            
            #pl.ax_steps.set_ylabel("$\\frac{d}{d\\epsilon} \\Gamma^{eii}$")

    pl.fig_steps.set_figheight(4.8)
    pl.fig_steps.set_figwidth(12.8)

    plt.tight_layout()
    plt.savefig(figure_output_dir + label +"_"+ str(slices[0]) + figures_ext)
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
        