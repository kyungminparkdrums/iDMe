import numpy as np
import awkward as ak

def cut0(events,info):
    name = "cut0"
    desc = "Dummy cut"
    plots = True
    return events, name, desc, plots

def cut1(events,info):
    name = "cut1"
    desc = "Nontrivial presel cuts"
    plots = True
    cut = (events.METFiltersFailBits == 0) & (events.trig.HLT_PFMET120_PFMHT120_IDTight == 1)
    if info['year'] == '2018':
        cut = cut & (~events.HEM.flag)
    events = events[cut]
    # require at least 1 good vtx
    events = events[ak.count(events.vtx.vxy,axis=1)>0]
    cut = ak.any((events.vtx.e1.pt>1) & (np.abs(events.vtx.e1.eta) < 2.4) & (events.vtx.e2.pt>1) & (np.abs(events.vtx.e2.eta) < 2.4),axis=1)
    events = events[cut]
    return events, name, desc, plots