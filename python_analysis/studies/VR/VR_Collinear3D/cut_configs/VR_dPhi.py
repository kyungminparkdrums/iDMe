import numpy as np
import awkward as ak
import sys
sys.path.append("../analysisTools/")
import analysisSubroutines as routines

def cut6(events,info):
    name = "cut6"
    desc = "Leading jet |eta| < 2.4"
    plots = False
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut7(events,info):
    name = "cut7"
    desc = "Leading jet pT > 80 GeV"
    plots = True
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "0.75 < dPhi(MET,leading jet) < 2.0"
    plots = True
    cut = (np.abs(events.PFJet.METdPhi[:,0]) < 2.0) & (np.abs(events.PFJet.METdPhi[:,0]) > 0.75)
    return events[cut], name, desc, plots