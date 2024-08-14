import numpy as np
import awkward as ak
from myHisto import myHisto

def make_histograms(info):
    h = myHisto()
    
    # matched reco electron, reco variables
    h.make('genPU_true','PU')
    h.make('genPU_obs','PU')

    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    if info["type"] == "signal" or info["type"] == "bkg":
        wgt = events.eventWgt/sum_wgt
    else:
        wgt = 1

    ### FILLING HISTOGRAMS ###
    
    #
    h.fill("genPU_true",pu=events.genPU.true,weight=wgt)
    h.fill("genPU_obs",pu=events.genPU.obs,weight=wgt)