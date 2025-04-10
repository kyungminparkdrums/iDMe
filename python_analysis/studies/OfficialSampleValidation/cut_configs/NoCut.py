import numpy as np
import awkward as ak
import sys
sys.path.append("../../../analysisTools/")
import analysisSubroutines as routines

def cut5(events,info):
    name = "cut0"
    desc = "Dummy"
    plots = True

    cut = (events.sel_vtx.sign == 1)|(events.sel_vtx.sign == -1)

    return events[cut], name, desc, plots