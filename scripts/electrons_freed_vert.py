# Comparing evolution of contributions by elements to freed electrons between multiple simulations.

import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': 8,
    'text.usetex': True,
    'pgf.rcfonts': False,
    'lines.linewidth':1.5,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook

DPI = 800

YLIM = [0,None]
#WIDTH = 3.49751
WIDTH = 3.484252*0.85 #7.24409436834/3 
FIGWIDTH = WIDTH
#FIGHEIGHT = WIDTH*3/4
FIGHEIGHT = max(WIDTH*3/4,2.2) # leave enough space for y axis label

XLIM = [None,None]

DASHED_AUGER = True

# ignore_CNO = False  # implemented in very hacky way
# CNO_only = False  # implemented in very hacky way



# ignore_CNO = False  # implemented in very hacky way
# CNO_only = True  # implemented in very hacky way

COMBINE_HEAVY_AND_LIGHT = False

SHOW_LEGEND = None# [False,True]
CUSTOM_LEGEND = None
#CUSTOM_LEGEND=("1s: CNO","2s: CNO","2p: CNO","1s: CNO,S,Gd","2s: CNO,S,Gd","2p: CNO,S,Gd")
#CUSTOM_LEGEND = ("Primary ionization only", "All ionization",)
def main():
    set_highlighted_excepthook()

    assert not (COMBINE_HEAVY_AND_LIGHT and (CNO_only == True or ignore_CNO == True))

    # Basic num arguments check
    assert len(sys.argv[1:]) > 0, "Usage: python compare_ion.py <sim_handle_1> <sim_handle_2> ..."
        
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
    label = data_folders[0]
    if len(data_folders) > 1:
         label+= "_"+data_folders[1]
    for ignore_CNO, CNO_only in [False,False],[False,True],[True,False]:
        for i, combined_primary in enumerate([False]):#enumerate([True,False]):
            
            ##The 3 graphs###
            if combined_primary and (ignore_CNO == True or CNO_only == True):
                continue
            if not combined_primary and (ignore_CNO == False and CNO_only == False):
                continue 

            
            ####
            if SHOW_LEGEND is not None:
                show_legend = SHOW_LEGEND[i]
            else:
                show_legend = True#combined_primary
            make_some_plots(data_folders,molecular_path,label,dname_Figures,combined_primary,distinguish_by_color_not_line=(not DASHED_AUGER),show_legend=show_legend,ignore_CNO=ignore_CNO, CNO_only=CNO_only)

def make_some_plots(mol_names,sim_output_parent_dir, label,figure_output_dir,combine_element_primary_electrons,distinguish_by_color_not_line=True,show_legend=True,ignore_CNO=False, CNO_only=False):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    assert not (CNO_only == True and ignore_CNO == True)  
    
    

    ############
    # File/directory names
    #######  
    figures_ext = ".png" #.png

       
    #cmap = plt.get_cmap("Dark2")
    cmap = plt.get_cmap("tab10")
    if distinguish_by_color_not_line:
        cmap = plt.get_cmap("tab20")
    split_continuums_to_load="all"
    e_cutoff = 500 # eV
    every = 1

    # Hacky way to get overlaid plots TODO
    pl = Plotter(mol_names[0],sim_output_parent_dir,split_continuums_to_load=split_continuums_to_load)
    
    pl.setup_axes(1)
    pl.fig.subplots_adjust(left=0.01,bottom=0.2,top=0.9,right=0.88)  # must not have tight bbox to avoid cutting off ylabel
    if ignore_CNO:
        pl.split_freeFiles = pl.split_freeFiles[6:] #hacky
    elif CNO_only:
        pl.split_freeFiles = pl.split_freeFiles[:6] #hacky
    continuums_set = []
    colors = []
    for c, continuum_fname in enumerate(pl.split_freeFiles):   

        if combine_element_primary_electrons:
            assert(len(pl.split_freeFiles)%2==0)
            if c%2==1:
                continue
            else:
                continuums_set.append([c,c+1])
        else:
            continuums_set.append([c])

        if ignore_CNO:
            c += 6 #hacky
        if distinguish_by_color_not_line:
            colors.append(cmap(c))  
        else:
            colors.append(cmap(int(c/2)))


    if COMBINE_HEAVY_AND_LIGHT:
        assert CNO_only is False
        assert ignore_CNO is False
        if combine_element_primary_electrons:
            new_continuums_set = [[],[]]
            for continuums in continuums_set[:3]:
                new_continuums_set[0].extend(continuums)
            for continuums in continuums_set[3:]:
                new_continuums_set[1].extend(continuums)
            

        else:
            new_continuums_set = [[],[],[],[]]
            for continuums in continuums_set[::2][:3]:
                new_continuums_set[0].append(continuums[0])
            for continuums in continuums_set[1::2][:3]:
                new_continuums_set[1].append(continuums[0])
            for continuums in continuums_set[::2][3:]:
                new_continuums_set[2].append(continuums[0])
            for continuums in continuums_set[1::2][3:]:
                new_continuums_set[3].append(continuums[0])
        continuums_set = new_continuums_set


    linestyles = ["solid"]*len(continuums_set)
    if not combine_element_primary_electrons and not distinguish_by_color_not_line:
        #linestyles = ["solid","dashed"]*int(len(continuums_set)/2)
        #param =(0,(2,1.5))
        param ="dashed"
        linestyles = ["solid",param]*int(len(continuums_set)/2)
        
    


    # if combined, each element of the continuums set has length of 2 (the continuums to be combined), otherwise length of 1.
    custom_legend = CUSTOM_LEGEND
    if custom_legend is None and COMBINE_HEAVY_AND_LIGHT:
        if combine_element_primary_electrons:
            custom_legend = ["Light atoms","Heavy atoms"]
        else: # NOTE this works because photo continuum files in Plotter class come before Auger, 
            custom_legend = ["Light atoms (Photo)","Light atoms (Auger)","Heavy atoms (Photo)","Heavy atoms (Auger)"]

        #custom_legend = ["Z<10","Z>=10"]

    def set_legend_wide_fig():
        ncol = 1
        #bbox_to_anchor = (1.02, 1)
        bbox_to_anchor = None
        if not combine_element_primary_electrons:
            ncol = 2
        if custom_legend is None: 
            ax.legend(bbox_to_anchor=bbox_to_anchor,loc='upper left', ncol=ncol,handlelength=1)  # Top right legend.
        else:
            handles,_ = ax.get_legend_handles_labels()
            handles = list(handles)
            assert len(custom_legend)==len(handles), f"{len(custom_legend)} != {len(handles)}"
            ax.legend(handles,custom_legend,bbox_to_anchor=bbox_to_anchor,loc='upper left', ncol=ncol,handlelength=1)  
    
    def set_horiz_legend():
        ncol = len(continuums_set)
        bbox_to_anchor = (-0.045,1.18)
        #bbox_to_anchor = None
        handles,labels = ax.get_legend_handles_labels()
        if custom_legend is not None: 
            handles = list(handles)
            assert len(custom_legend)==len(handles), f"{len(custom_legend)} != {len(handles)}"
            labels = custom_legend
        ax.legend(handles,labels,loc='upper left',bbox_to_anchor=bbox_to_anchor, ncol=ncol,borderpad=0.2,handletextpad=0.2,handlelength=0.3,columnspacing=0.35) 
    
    #2col legend
    def set_phot_aug_legend(photo_handles,aug_handles,loc= "upper left",): #sep,hstart=0,vstart=1,
        actual_handles,actual_labs = ax.get_legend_handles_labels()
        elements = actual_labs[::2]
        for i, e in enumerate(elements):
            elements[i] = e.split(' ')[0]

        assert(len(photo_handles) == len(aug_handles))
        i = 0
        legends = []
        # for photo, aug,element in zip (photo_handles,aug_handles,elements):
        #     process_labels=[None,None]
        #     if i == 0:
        #         process_labels = ["Photo.","Auger"]
        #     print(process_labels,element)
        #     leg = plt.legend(handles=(photo,aug),labels=[process_labels],title=element,loc = loc,bbox_to_anchor=[hstart+sep*i, vstart],
        #             ncol=1,
        #             borderpad=0.2,handletextpad=-0.5,handlelength=0.3,columnspacing=0.35
        #                     )
            
        #     legends.append(leg)
        #     i+=1

        # gap_a = 5
        # labels = ['Photo.' + '\\ '*gap_a,'Auger'+ '\\ '*gap_a]
        # labels.extend([None]*(len(actual_handles)-2))
        # gap_title = 7
        # gap_title_left = 15
        # element_title = '\\ '*gap_title_left
        # for element in elements:
        #     element_title += element + '\\ '*gap_title
        # print(element_title)
        # leg = plt.legend(actual_handles,labels,markerfirst = False, ncol=int(len(actual_handles)/2),loc=loc,
        #                  title=element_title,
        #                  borderpad=0.2,handletextpad=0,handlelength=1.5,columnspacing=1,
        #                  )
        gap_a = 2
        labels = []
        for element in elements:
            labels.append(element)
        labels.extend([None]*len(elements))
        gap_title = 1

        title = 'Phot.'+'\\ '*gap_title+'Aug.'
        print(title)
        handles = []
        handles.extend(photo_handles)
        handles.extend(aug_handles)
        leg = plt.legend(handles,labels,markerfirst = False, ncol=2,loc=loc,
                         title=title,
                         borderpad=0.3,handletextpad=1,handlelength=1.5,columnspacing=0.65,
                         bbox_to_anchor=(-0.02,1.03),
                         framealpha=0,
                         )
        leg._legend_box.align = "right"
        leg.get_frame().set_linewidth(0)
        #leg.get_frame().set_facecolor('#e6e6e6') #'whitesmoke' 'snow'
        #leg.get_frame().set_facecolor('#e6e6e6') #'whitesmoke' 'snow'
        #leg.get_frame().set_facecolor('#e6e6e6') #'whitesmoke' 'snow'
        #leg.get_frame().set_boxstyle('Square')
        legends.append(leg)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        for l in legends:
            plt.gca().add_artist(l)
        #table = ax.table()
        # leg = plt.legend(line_columns, ['']*len(line_columns), 
        #      title='No Prop    Prop +    Prop *',  
        #      ncol=3, numpoints=1, handletextpad=-0.5)

    handles = []
    for i, continuums in enumerate(continuums_set):    
        pl.update_free
        pl.num_plotted = 0 # ƪ（˘へ˘ ƪ）
        show_pulse_profile = True
        if i != 0:
            show_pulse_profile = False
        _, subset_of_handles = pl.plot_electrons_freed(continuums,e_cutoff=e_cutoff,every=every,ylim=YLIM,show_pulse_profile=show_pulse_profile,color=colors[i],ls=linestyles[i],
                                pulse_profile_height_factor=0.87)
        ax = pl.axs[0][0]
        handles.extend(subset_of_handles)

    if show_legend:
        if combine_element_primary_electrons:
            set_horiz_legend()
            #set_legend_wide_fig() 
        else:
            set_phot_aug_legend(handles[::2],handles[1::2])   

    plt.gcf().set_figwidth(FIGWIDTH)
    plt.gcf().set_figheight(FIGHEIGHT)
    #plt.gcf().tight_layout()
    #plt.tight_layout()




    qualifier = f"_electrons_freed"
    if combine_element_primary_electrons:
        qualifier += "_comb"
    if DASHED_AUGER and not combine_element_primary_electrons:
        qualifier+="_dashed_aug" 
    if ignore_CNO:
        qualifier += "_no_CNO"
    if CNO_only:
        qualifier += "_CNO_only"


    if COMBINE_HEAVY_AND_LIGHT:
        qualifier += "_heavy_and_light"

    #bbox = pl.fig.gca().get_tightbbox(for_layout_only=False)
    #plt.gca().set_axis_off()
    #plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
                #hspace = 0, wspace = 0)
    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    #plt.margins(0,0)


    plt.savefig(figure_output_dir + label +qualifier + figures_ext, dpi=DPI)


    #plt.savefig(figure_output_dir + label +qualifier + figures_ext,bbox_inches='tight',dpi=DPI)

    plt.close()
    print("Done")

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option