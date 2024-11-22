import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Preselection"
    plots = False
    return events, name, desc, plots

def cut1(events,info):
    name = "cut1"
    desc = "Njets < 3"
    plots = False
    cut = ak.count(events.PFJet.pt,axis=1) < 3
    return events[cut], name, desc, plots

def cut2(events,info):
    name = "cut2"
    desc = "Leading jet |eta| < 2.4"
    plots = False
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "Leading jet pT > 80 GeV"
    plots = False
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "dPhi(MET,leading jet) > 2.0"
    plots = False
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 2.0
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = True
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots
