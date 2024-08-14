import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Dummy cut"
    plots = True
    return events, name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "dPhi(MET,leading jet) < 1.5"
    plots = True
    cut = np.abs(events.PFJet.METdPhi[:,0]) < 1.5
    return events[cut], name, desc, plots