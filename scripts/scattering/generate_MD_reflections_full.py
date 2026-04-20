#%%
import sys
sys.path.append('/home/speno/AC4DC')

from generate_MD_reflections import generate_MD_reflections

#TODO Make importable: convert from command line script.



SPI=False


plasma_sim_handle=sys.argv[1]
MD_results_parent_dir=sys.argv[2]
timespan_ps_MD=float(sys.argv[3])

#ground_truth_pdb=None if len(sys.argv)<4 else sys.argv[3]
#assert ground_truth_pdb is not None # XXX
ground_truth_pdb=sys.argv[4]


start_time_ps = 0 if len(sys.argv)<6 else float(sys.argv[5])

extra_tag = "" if len(sys.argv)<7 else sys.argv[6]

base_tag = "All"
electronic_damage = True
nuclear_damage = True
generate_MD_reflections(ground_truth_pdb,plasma_sim_handle,MD_results_parent_dir,timespan_ps_MD,base_tag,
extra_tag=extra_tag,start_time_ps=start_time_ps,
electronic_damage=electronic_damage,nuclear_damage=nuclear_damage,SPI=SPI)
#all_reflections_to_scalepack() 




# %%
