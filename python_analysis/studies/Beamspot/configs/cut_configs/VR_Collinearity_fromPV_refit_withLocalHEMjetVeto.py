import numpy as np
import awkward as ak
import sys

'''
def cut1(events,info):
    name = "cut1"
    desc = "Preselections (no HEM veto) & n(good_vtx) & n(Jet) Cuts"
    plots = True

    return events, name, desc, plots
'''

def cut2(events,info):
    name = "cut2"
    desc = "HEM Jet Veto (eta upper bound bug)"
    plots = True

    if len(events) != 0:
        cut = events.hasHEMjetBug == 0
        #print(f'HEM Jet Veto Pass (with eta bug): {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []

    return events[cut], name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "HEM Jet Veto (additionally veto missing eta region)"
    plots = True

    if len(events) != 0:
        cut = events.hasHEMjet == 0
        #print(f'HEM Jet Veto Pass (additional eta region): {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []

    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "HEM electron Veto"
    plots = True

    if len(events) != 0:
        cut = (events.hasHEMelecPF == 0) & (events.hasHEMelecLpt == 0)
        #print(f'HEM electron Veto Pass: {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []

    return events[cut], name, desc, plots

def cut5(events,info):
    # UL b-tag threshold recommendations from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    # using the medium WP, as in Andre's version of iDM
    name = "cut5"
    desc = "No b-tagged jets"
    plots = True
    bTag = events.PFJet.bTag
    # DeepFlavour working points for UL samples
    if info["year"] == 2018:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        wp = 0.0490 # loose
        #wp = 0.2783 # medium
        #wp = 0.7100 # tight
    if info["year"] == 2017:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
        wp = 0.0532 # loose
        #wp = 0.3040 # medium
        #wp = 0.7476 # tight
    if info["year"] == "2016_preVFP":
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP
        wp = 0.0508 # loose
        #wp = 0.2598 # medium
        #wp = 0.6502 # tight
    if info["year"] == "2016_postVFP":
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP
        wp = 0.0480 # loose
        #wp = 0.2489 # medium
        #wp = 0.6377 # tight
    pass_bTag = bTag > wp
    n_bTag_Jets = ak.count(events.PFJet[pass_bTag].pt,axis=1)
    cut = n_bTag_Jets == 0
    return events[cut], name, desc, plots

def cut6(events,info):
    name = "cut6"
    desc = "Leading jet |eta| < 2.4"
    plots = True
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
    desc = "dPhi(MET,leading jet) > 1.5"
    plots = True
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 1.5
    return events[cut], name, desc, plots

def cut9(events,info):
    name = "cut9"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = True
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots

def cut10(events,info):
    name = "cut10"
    desc = "Vtx cos_collinear_fromPV_refit < 0"
    plots = True
    cut = events.sel_vtx.cos_collinear_fromPV_refit < 0
    
    return events[cut], name, desc, plots

'''
def cut11(events,info):
    name = "cut11"
    desc = "Vtx cos_collinear_fromPV_refit > -0.99"
    plots = True

    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.99

    return events[cut], name, desc, plots

def cut12(events,info):
    name = "cut12"
    desc = "Vtx cos_collinear_fromPV_refit > -0.98"
    plots = False

    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.98

    return events[cut], name, desc, plots

def cut13(events,info):
    name = "cut13"
    desc = "Vtx cos_collinear_fromPV_refit > -0.97"
    plots = False

    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.97

    return events[cut], name, desc, plots

def cut14(events,info):
    name = "cut14"
    desc = "Vtx cos_collinear_fromPV_refit > -0.96"
    plots = False

    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.96

    return events[cut], name, desc, plots

'''
def cut11(events,info):
    name = "cut11"
    desc = "Vtx cos_collinear_fromPV_refit > -0.95"
    plots = True
    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.95

    mask_cut_fail = ~cut
    print('Cos(collinear) [-1, -0.95] event #: ', events.eventNum[mask_cut_fail])
    
    return events[cut], name, desc, plots

'''
def cut16(events,info):
    name = "cut16"
    desc = "Vtx cos_collinear_fromPV_refit > -0.9"
    plots = False
    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.9
    
    return events[cut], name, desc, plots

def cut17(events,info):
    name = "cut17"
    desc = "Vtx cos_collinear_fromPV_refit > -0.85"
    plots = False
    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.85
    
    return events[cut], name, desc, plots

def cut18(events,info):
    name = "cut18"
    desc = "Vtx cos_collinear_fromPV_refit > -0.8"
    plots = False
    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.8
    
    return events[cut], name, desc, plots

def cut19(events,info):
    name = "cut19"
    desc = "Vtx cos_collinear_fromPV_refit > -0.75"
    plots = False
    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.75
    
    return events[cut], name, desc, plots

def cut20(events,info):
    name = "cut20"
    desc = "Vtx cos_collinear_fromPV_refit > -0.7"
    plots = False
    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.7
    
    return events[cut], name, desc, plots


def cut11(events,info):
    name = "cut11"
    desc = "Lxy (from PV) > 0.1"
    plots = True

    vx_fromPV = events.sel_vtx.vx - events.PV.x
    vy_fromPV = events.sel_vtx.vy - events.PV.y

    vxy_fromPV = np.sqrt(vx_fromPV * vx_fromPV + vy_fromPV * vy_fromPV)
    
    cut = vxy_fromPV > 0.1
    
    return events[cut], name, desc, plots
'''

'''
def cut11(events,info):
    name = "cut11"
    desc = "BDT"
    plots = True

    thres = 0.96 # BDT threshold; inference is done in analysisSubroutines

    if len(events) != 0:
        cut = events['BDTScore'] > thres
    else:
        cut = [] # this if-else is a workaround of some weird bug
    
    return events[cut], name, desc, plots
'''