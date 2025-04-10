import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    # electron
    h.make("PFJet_bTagEffNumPt",'ele_pt')
    h.make("PFJet_bTagEffNumPt",'ele_pt')
    
    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    if info["type"] == "signal" or info["type"] == "bkg":
        wgt = events.eventWgt/sum_wgt
    else:
        wgt = 1

    sel_vtx = events.sel_vtx
    ### FILLING HISTOGRAMS ###
    
    h.fill("sel_e1_pt",pt=sel_vtx.e1.pt,weight=wgt)
    