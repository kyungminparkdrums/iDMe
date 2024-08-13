import numpy as np
import awkward as ak
import sys

'''
def cut0(events,info):
    name = "cut0"
    desc = "No Cut"
    plots = True
    return events, name, desc, plots
'''

def cut1(events,info):
    name = "cut1"
    desc = "Pass MET Filters"
    plots = False
    cut = events.METFiltersFailBits == 0
    return events[cut], name, desc, plots

def cut2(events,info):
    name = "cut2"
    desc = "HEM Veto"
    plots = False
    if info["year"] == 2018:
        cut = ~events.HEM.flag
        return events[cut], name, desc, plots
    else:
        return events, name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "MET Trigger (120 GeV)"
    plots = True

    '''
    if info["year"] == 2016:
        cut = (events.fired16 & (1<<9)) == (1<<9)
    if info["year"] == 2017:
        cut = (events.fired17 & (1<<9)) == (1<<9)
    if info["year"] == 2018:
        cut = (events.fired18 & (1<<13)) == (1<<13)
    '''
    cut = events.trig.HLT_PFMET120_PFMHT120_IDTight == 1
    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "MET > 200 GeV"
    plots = True
    cut = events.PFMET.pt > 200
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "0 < nJets < 3 (pT > 30 GeV)"
    plots = True
    nJets = ak.count(events.PFJet.pt,axis=1)
    cut = (nJets > 0) & (nJets < 3)
    return events[cut], name, desc, plots

def cut6(events,info):
    # UL b-tag threshold recommendations from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    # using the medium WP, as in Andre's version of iDM
    name = "cut6"
    desc = "No b-tagged jets"
    plots = False
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