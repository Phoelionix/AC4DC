#%%
import sys
sys.path.append('/home/speno/AC4DC')

from generate_MD_reflections import generate_MD_reflections

#TODO Make importable: convert from command line script.



SPI=False


plasma_sim_handle=sys.argv[1]
MD_results_parent_dir=sys.argv[2]

#ground_truth_pdb=None if len(sys.argv)<4 else sys.argv[3]
#assert ground_truth_pdb is not None # XXX
ground_truth_pdb=sys.argv[3]

extra_tag = "" if len(sys.argv)<5 else sys.argv[4]

base_tag = "All"
electronic_damage = True
nuclear_damage = True
generate_MD_reflections(ground_truth_pdb,plasma_sim_handle,MD_results_parent_dir,base_tag,
extra_tag=extra_tag,electronic_damage=electronic_damage,nuclear_damage=nuclear_damage,SPI=SPI)
#all_reflections_to_scalepack() 




# %%
