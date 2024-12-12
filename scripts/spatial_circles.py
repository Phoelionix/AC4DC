from spatial_interactive_core import SpatialInteractive
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook




target_handles = sys.argv[1:]  
if  len(target_handles) != 1:
    print("Usage: python3 " +path.basename(__file__)+" Carbon_1")
    exit()

molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
dname_Figures = "../../output/_Graphs/interactives/"
dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"

######

SI = SpatialInteractive(target_handles,molecular_path)
SI.plot_circles("C",1,6)
SI.save_interactive(dname_Figures,"circle_test")
