import numpy as np
import awkward as ak
import sys
import json
import pandas as pd
import os

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
    desc = "|dPhi(MET,leading jet)| > 1.5"
    plots = True
    cut = np.abs(events.PFJet.METdPhi[:,0]) > 1.5
    return events[cut], name, desc, plots

def cut9(events,info):
    name = "cut9"
    desc = "|dPhi(MET,all jets)| > 0.75"
    plots = True
    cut = ak.all(np.abs(events.PFJet.METdPhi) > 0.75,axis=1)
    return events[cut], name, desc, plots

def cut10(events,info):
    name = "cut10"
    desc = "Vtx cos(collinear) < 0"
    plots = True
    cut = events.sel_vtx.cos_collinear_fromPV_refit < 0
    
    return events[cut], name, desc, plots

def cut11(events,info):
    name = "cut11"
    desc = "Vtx cos(collinear) > -0.99"
    plots = True
    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.99

    #mask_cut_fail = ~cut

    #evt_dict = {}
    #evt_dict['dataset'] = events.metadata["dataset"]
    #evt_dict['runNum'] = events.runNum[mask_cut_fail].to_list()
    #evt_dict['LS'] = events.lumiSec[mask_cut_fail].to_list()
    #evt_dict['eventNum'] = events.eventNum[mask_cut_fail].to_list()
    
    #print(evt_dict)

    #if 'data' in evt_dict['dataset']:
    #    df = pd.DataFrame(evt_dict)
    #    output_path = 'cos_minus_0p99_1_events_data.csv'
    #    df.to_csv(output_path, mode='a', header=not os.path.exists(output_path))
    #else:
    #    df = pd.DataFrame(evt_dict)
    #    output_path = 'cos_minus_0p99_1_events_bkg.csv'
    #    df.to_csv(output_path, mode='a', header=not os.path.exists(output_path))

    
    return events[cut], name, desc, plots

'''
def cut11(events,info):
    name = "cut11"
    desc = "dPhi(vtx, MET)"
    #desc = "Vtx log10(dxy/dz) > -1.25"
    plots = True
    #cut = (np.log10(np.abs(events.sel_vtx.e1.dxy/events.sel_vtx.e1.dz)) > -1.25) & (np.log10(np.abs(events.sel_vtx.e2.dxy/events.sel_vtx.e2.dz)) > 1.25)
  
    return events[cut], name, desc, plots
'''

'''
def cut12(events,info):
    name = "cut12"
    desc = "Vtx cos(collinear) > -0.9"
    plots = True
    cut = events.sel_vtx.cos_collinear_fromPV_refit > -0.9

    mask_cut_fail = ~cut

    evt_dict = {}
    evt_dict['dataset'] = events.metadata["dataset"]
    evt_dict['runNum'] = events.runNum[mask_cut_fail].to_list()
    evt_dict['LS'] = events.lumiSec[mask_cut_fail].to_list()
    evt_dict['eventNum'] = events.eventNum[mask_cut_fail].to_list()
    
    print(evt_dict)

    df = pd.DataFrame(evt_dict)
    output_path = 'cos_minus1_events.csv'
    df.to_csv(output_path, mode='a', header=not os.path.exists(output_path))
    
    return events[cut], name, desc, plots
'''