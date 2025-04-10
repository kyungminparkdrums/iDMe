import numpy as np
import awkward as ak
import sys
sys.path.append("../../../analysisTools/")
import analysisSubroutines as routines
import vector

def cut1(events,info):
    name = "cut1"
    desc = "Dummy"
    plots = True
    cut = events.METFiltersFailBits == 0
    
    return events[cut], name, desc, plots
