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

def cut2(events,info):
    name = "cut2"
    desc = "Signal reconstructed"
    plots = True
    cut = events.signalReconstructed == 1
    
    return events[cut], name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "Signal reco LL"
    plots = True

    reco_type = ak.to_numpy(np.stack((events.GenEle.matchType,events.GenPos.matchType),axis=1))
    
    mask_reco_R = reco_type == 'R'
    mask_reco_L = reco_type == 'L'
    mask_reco_None = reco_type == 'None'

    mask_reco_RR = ak.values_astype(ak.sum(mask_reco_R, axis=1) == 2, int)
    mask_reco_LL = ak.values_astype(ak.sum(mask_reco_L, axis=1) == 2, int)
    mask_reco_None = ak.values_astype(ak.sum(mask_reco_None, axis=1) > 0, int)
    mask_reco_RL = ak.values_astype((ak.sum(mask_reco_R, axis=1) == 1)&(ak.sum(mask_reco_L, axis=1) == 1), int)

    cut = mask_reco_LL == 1

    return events[cut], name, desc, plots