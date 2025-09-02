#%%
import matplotlib
#matplotlib.use("pgf")#
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'axes.titlesize':8,     # fontsize of the axes title
    'axes.labelsize':8,    # fontsize of the x and y labels   
    'ytick.labelsize':8,
    'xtick.labelsize':8,
    'legend.fontsize':6,    
    'legend.title_fontsize':8,    
    'lines.linewidth':1,
})
from scatter import *
from multiprocessing import Pool
from sample_scalepack import all_reflections_to_scalepack
import inspect
import datetime
from core_functions import get_sim_params

#TODO auto generate ideal, undamaged, damaged. not undamaged and damaged. (undamaged is called ideal)

NUM_PARALLEL=1
SEEDED = False
RANDOM_WATER = False; NUM_RANDOM_WATER = 0 # 702
DEBUG_WATER = False
QUICK_TEST = False
SKIP_UNDAMAGED = False


target_handle =  "copper_sulfate_above_e12_15" #"copper_sulfate_below_e13_3#"copper_sulfate_above_e12_long_1"#"copper_sulfate_above_e12_14"
#target_handle = "copper_sulfate_above_e12_15" #"copper_sulfate_below_e13_3#"copper_sulfate_above_e12_long_1"#"copper_sulfate_above_e12_14"

# Ideal intensities purely as an extremely inefficient lazy way to scale by pulse profile
def plot_spots_at_time_slices(Intensities, ideal_intensities, miller_indices,times,sqrt=False,scale_by_pulse_profile=True):
    plt.close()
    stringified = []
    Intensities = copy.deepcopy(Intensities)
    scale_I = np.sum(Intensities[0])
    for elem in miller_indices:
        stringified.append(str(elem[0])+str(elem[1])+str(elem[2])) 
    for I,I_ideal,col, w,  in zip(Intensities,ideal_intensities,['b','g','r'],[1,0.7,0.4]):     
        if scale_by_pulse_profile:
            I *= scale_I/np.sum(I_ideal) #normalise
        plt.bar(stringified,np.sqrt(I) if sqrt else I,alpha=1,color=col,width=w)
    plt.ylim(0,np.sqrt(np.max(Intensities)))
    plt.xticks(rotation="vertical")
    #plt.show()

def plot_biggest_spots_over_time(Intensities, ideal_intensities, miller_indices,times,xlim=[None,None],ylim=[None,None],sqrt=False,scale_by_pulse_profile=False,num_points=4,save_path = None,title="",true_normed=False,structure_factors=False):
    plt.close()
    Intensities = copy.deepcopy(Intensities)
    
  
    scale_I = np.sum(ideal_intensities[0])
    # initialise dict of reflection intensities 
    data = {}
    hkl_keys=[] #ordered keys
    for elem in miller_indices:
        k=str(elem[0])+str(elem[1])+str(elem[2])
        data[k]=[]
        hkl_keys.append(k)
    # square=True
    # if square:
    #     times = np.insert(times,0,0)
    #     Intensities = np.insert(Intensities,0,ideal_intensities[i])
    for i in range(len(Intensities)): # for each time step
        I = Intensities[i]; I_ideal = ideal_intensities[i]
        assert(np.sum(I)>0)
        if scale_by_pulse_profile:
            I *= scale_I/np.sum(I_ideal)
        if true_normed:
            I/=np.sum(I)
        for hkl, v in zip(hkl_keys,I):
            data[hkl].append(np.sqrt(v) if sqrt else v)
    
        #plt.bar(stringified,np.sqrt(I),alpha=1,color=col,width=w)
    fig_width = 3.49751*0.75 # 20
    fig_height = fig_width*3/4 # 20
    fig,ax = plt.subplots(figsize=(fig_width,fig_height))
    

    Y_list =np.array(list(data.values()))
    if structure_factors:
        Y_list = np.sqrt(Y_list)
    sorts = np.argsort(Y_list[:,0])[::-1]
    hkl_list = np.array(list(data.keys()))[sorts][num_points[0]:num_points[1]]
    Y_list = Y_list[sorts][num_points[0]:num_points[1]]

    diff = (Y_list[:,0]-Y_list[:,-1])/Y_list[:,0]
    print("max diff",np.max(diff))
    print("min diff",np.min(diff))


    for hkl, Y in zip(hkl_list,Y_list):
        #ax.scatter(times,Y,label=hkl,s=4)
        #ax.plot(times,Y,'-o',label=hkl,markersize=2.5)
        ax.plot(times,Y,'-|',label=hkl,markersize=2)
        #A = np.vstack([times, np.ones(len(times))]).T
        #m,c =  np.linalg.lstsq()
        
    #ax.set_yscale("log")
    leg =ax.legend(loc="upper left",bbox_to_anchor=[1,1],ncol=2,title="$hkl\\:$",
        borderpad=0.2,handletextpad=0.2,handlelength=0.5,columnspacing=0.35) 
   
    leg._legend_box.align = "right"

    ax.set_ylabel("$$I_{hkl}$$")
    if structure_factors:
        ax.set_ylabel("$$\\left|F(hkl)\\right|$$")
    if scale_by_pulse_profile:
        #ax.set_ylabel("$$\\frac{I}{\\Omega}$$")
        ax.set_ylabel("$$I_{hkl}/\\Omega$$")
        if structure_factors:
            ax.set_ylabel("$$PlaceH$$")
    if true_normed:
        #ax.set_ylabel("$$\\frac{I}{\\Omega}$$")
        ax.set_ylabel("$$I_{hkl}/\\Sigma I$$")
        if structure_factors:
            ax.set_ylabel("$$PlaceH$$")
    assert(not sqrt)
    ax.set_xlabel("Time (fs)")
    ax.set_title(title)
    #ax.set_xticks((0,1,2,3,4,5))
    fig.subplots_adjust(left=0.15, bottom=0.1, right=0.8, top=0.93, wspace=0.3, hspace=0.4) # 3 orb density plots row

    ax.set_ylim(*ylim)
    ax.set_xlim(*xlim)
    #plt.ylim(0,None)

    plt.tight_layout(pad=0)

    if save_path is not None:
        plt.savefig(save_path,dpi=300)
    #plt.show()




if __name__ == "__main__":
    import glob
    #target_handle = "copper_sulfate_above_e12_15"
    src_file_path = inspect.getfile(lambda: None)
    scattering_dir = path.abspath(path.join(src_file_path ,"../"))+"/"
    RESULTS_LOCAL_PATH = "results/"


    #num_results = len(times)*(1+(SKIP_UNDAMAGED==False))
    #OldestToLatest = sorted(glob.glob(os.path.join(scattering_dir+RESULTS_LOCAL_PATH, '*/')), key=os.path.getmtime)
    dir = scattering_dir+RESULTS_LOCAL_PATH+target_handle+"_real"
    #print (glob.glob(os.path.join(dir,"*/")))
    #    #for file in glob.glob(dir+"/*"):
    Intensities = []
    ideal_intensities=[]
    files = []
    subdirs =[]
    for filepath in os.listdir(dir):
        subdir = dir+"/"+filepath
        for file in os.listdir(subdir):
            for i, c in  enumerate(filepath):
                if c== "-":
                    break
            files.append(file)
            subdirs.append(subdir)
    assert len(files)==1
    
    for file,subdir in zip(files,subdirs):
        compare_dir=subdir.replace("_real","_ideal")
        dmged_result,undmged_result = get_result(file,subdir,compare_dir=compare_dir)
        Intensities = dmged_result.I[:]
        #Intensities = dmged_result.I[1:-1]
        #ideal_intensities = undmged_result.I[1:-1]
        times = dmged_result.T[:]
        #times = dmged_result.T[1:-1]
        miller_indices = dmged_result.miller_indices
    #plot_spots_at_time_slices(Intensities,ideal_intensities,miller_indices,times)
    dir_figures = "../../../output/_Graphs/plots/"
    dir_figures = path.abspath(path.join(__file__ ,dir_figures)) + "/"

    param_dict, param_name_list,unit_list,pulse_profile = get_sim_params(
        target_handle,get_intensities_at_times=np.array(times))
    energy,fwhm,photon_count = param_dict["energy"],param_dict["width"],param_dict["fluence"]

    #note = f"Decay of strongest {num_points} Bragg peaks"
    #title = f"Pulse: "+"$10_power_photon_coun$"+f" ${energy/1000}$"+" keV ph $\\cdot$µm$^"+"{-2}$"+f", ${fwhm}$ fs FWHM"
    title = f"Pulse: "+"$10^{12}$"+f" ${energy/1000}$"+" keV ph $\\cdot$µm$^"+"{-2}$"+f", ${fwhm}$ fs FWHM"
    times+=5
   # plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,ylim=[2e4,4e4],sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"_high_rfl.png",title=title)
   # plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,ylim=[0,1.2e4],sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"_low_rfl.png",title=title)
    xlim=[0,5]
    plot_biggest_spots_over_time(Intensities,pulse_profile,miller_indices,times,ylim=[0,None],sqrt=False,
                                 num_points=[0,10],xlim=xlim,save_path=dir_figures+target_handle+"_1rfl.png",title=title)
    plot_biggest_spots_over_time(Intensities,pulse_profile,miller_indices,times,ylim=[0.8,None],sqrt=False,
                                 num_points=[10,20],xlim=xlim,save_path=dir_figures+target_handle+"_2rfl.png",title=title)
    plot_biggest_spots_over_time(Intensities,pulse_profile,miller_indices,times,ylim=[0,None],sqrt=False,
                                 num_points=[20,30],xlim=xlim,save_path=dir_figures+target_handle+"_3rfl.png",title=title)
    plot_biggest_spots_over_time(Intensities,pulse_profile,miller_indices,times,ylim=[0,None],sqrt=False,
                                 num_points=[0,10],xlim=xlim,save_path=dir_figures+target_handle+"_1rfl_normed.png",title=title,
                                 true_normed=True)
    plot_biggest_spots_over_time(Intensities,pulse_profile,miller_indices,times,ylim=[0,None],sqrt=False,
                                 num_points=[0,10],xlim=xlim,save_path=dir_figures+target_handle+"_1rfl_scaled_by_profile.png",title=title,
                                 scale_by_pulse_profile=True)

    #plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,scale_by_pulse_profile=True,sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"pulse_scaled_rfl.png",title=title)
    #plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,scale_by_pulse_profile=True,sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"normed_rfl.png",title=title,true_normed=True)
    assert len(files)>0
    
    # for n in range(num_results):
    #     all_reflections_to_scalepack(
    #         OldestToLatest[-n-1].split("/")[-2],
    #         scattering_dir+RESULTS_LOCAL_PATH,
    #         out_dir="/home/speno/PhenixWorkspace/data/",
    #         tag_override=""
    #     )  


# %%
from core_functions import get_sim_params



# %%
