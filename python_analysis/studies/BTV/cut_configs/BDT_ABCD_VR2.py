import numpy as np
import awkward as ak
import sys
sys.path.append("../../../analysisTools/")
import analysisSubroutines as routines


def cut5(events,info):
    # UL b-tag threshold recommendations from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    # using the medium WP, as in Andre's version of iDM
    name = "cut5"
    desc = "No b-tagged jets"
    plots = False
    bTag = events.PFJet.bTag
    # DeepFlavour working points for UL samples
    if info["year"] == 2018:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        #wp = 0.0490 # loose
        wp = 0.2783 # medium
        #wp = 0.7100 # tight
    if info["year"] == 2017:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
        #wp = 0.0532 # loose
        wp = 0.3040 # medium
        #wp = 0.7476 # tight
    if info["year"] == "2016_preVFP":
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP
        #wp = 0.0508 # loose
        wp = 0.2598 # medium
        #wp = 0.6502 # tight
    #if info["year"] == "2016_postVFP":
    if info["year"] == 2016:
        # twiki : https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP
        #wp = 0.0480 # loose
        wp = 0.2489 # medium
        #wp = 0.6377 # tight
    pass_bTag = bTag > wp
    n_bTag_Jets = ak.count(events.PFJet[pass_bTag].pt,axis=1)
    cut = n_bTag_Jets == 0
    return events[cut], name, desc, plots


def cut6(events,info):
    name = "cut6"
    desc = "Leading jet pT > 80 GeV"
    plots = True
    cut = events.PFJet.pt[:,0] > 80
    return events[cut], name, desc, plots


def cut7(events,info):
    name = "cut7"
    desc = "Leading jet |eta| < 2.4"
    plots = False
    cut = np.abs(events.PFJet.eta[:,0]) < 2.4
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "dPhi(MET,leading jet) > 2.0"
    plots = False
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 2.0
    return events[cut], name, desc, plots

def cut9(events,info):
    name = "cut9"
    desc = "dPhi(MET,all jets) > 0.75"
    plots = True
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots

def cut10(events,info):
    name = "cut10"
    desc = "OSSF"
    plots = True
    cut = events.sel_vtx.sign == -1
    return events[cut], name, desc, plots

'''
def cut11(events,info):
    name = "cut11"
    desc = "theta < 90"
    plots = True
    cut = np.degrees(np.arccos(events.sel_vtx.cos_collinear_fromPV_refit)) < 90
    return events[cut], name, desc, plots
'''

def cut11(events,info):
    name = "cut11"
    desc = "BDT"
    plots = True

    #thres = 0.94
    thres_low = 0.2
    #thres_upp = 0.45
    thres_upp = 0.5
    
    if len(events) != 0:
        input = routines.makeBDTinputs_ABCD(events)
    
        model = './models/BDT_NJetG0_allSigns_goodvtx_v11_max_depth_7_n_estimators_800_lr_0.005.json'
        score_BDT = routines.getBDTscore(input, model)

        #cut = score_BDT > thres
        cut = (score_BDT < thres_upp) & (score_BDT > thres_low) 
        
        print(f'BDT Pass: {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []
    
    return events[cut], name, desc, plots
