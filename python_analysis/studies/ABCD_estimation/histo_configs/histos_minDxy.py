import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    h.make("sel_e1_dxy",'dxy')
    h.make("sel_e2_dxy",'dxy')

    h.make("sel_min_dxy",'dxy')
 
    h.make("sel_e1_dxy_refit", 'dxy')
    h.make("sel_e2_dxy_refit", 'dxy')

    h.make("sel_min_dxy_refit",'dxy')
    
    return h

subroutines = []

def dxy_custom(vx, vy, beamspot_x, beamspot_y, px, py, pt):
    dxy = ((-(vx - beamspot_x) * py) + (vy - beamspot_y)*px)/pt

    return dxy

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    if info["type"] == "signal" or info["type"] == "bkg":
        wgt = events.eventWgt/sum_wgt
    else:
        wgt = 1

    sel_vtx = events.sel_vtx
    ### FILLING HISTOGRAMS ###

    h.fill("sel_e1_dxy",dxy=np.abs(sel_vtx.e1.dxy),weight=wgt)
    h.fill("sel_e2_dxy",dxy=np.abs(sel_vtx.e2.dxy),weight=wgt)

    h.fill("sel_min_dxy",dxy=np.minimum(np.abs(sel_vtx.e1.dxy),np.abs(sel_vtx.e2.dxy)),weight=wgt)

    h.fill("sel_e1_dxy_refit", dxy=np.abs(sel_vtx.e1.refit_dxy), weight=wgt)
    h.fill("sel_e2_dxy_refit", dxy=np.abs(sel_vtx.e2.refit_dxy), weight=wgt)
    
    h.fill("sel_min_dxy_refit",dxy=np.minimum(np.abs(sel_vtx.e1.refit_dxy),np.abs(sel_vtx.e2.refit_dxy)),weight=wgt)