import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    #h.make("bdtscore_vs_sel_vtx_thetaColl_fromPV_refit",'bdtscore_zoom','theta')
    h.make("bdtscore_vs_sel_vtx_thetaColl_fromPV_refit_rad",'bdtscore_zoom','theta_rad')
    #h.make("bdtscore_vs_sel_vtx_thetaColl_fromPV",'bdtscore_zoom','theta')

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

    #h.fill("bdtscore_vs_sel_vtx_thetaColl_fromPV_refit",bdtscore=events.BDTScore,\
    #       theta=np.degrees(np.arccos(sel_vtx.cos_collinear_fromPV_refit)),weight=wgt)
    h.fill("bdtscore_vs_sel_vtx_thetaColl_fromPV_refit_rad",bdtscore=events.BDTScore,\
           theta_rad=np.arccos(sel_vtx.cos_collinear_fromPV_refit),weight=wgt)
    #h.fill("bdtscore_vs_sel_vtx_thetaColl_fromPV",bdtscore=events.BDTScore,\
    #       theta=np.degrees(np.arccos(sel_vtx.cos_collinear_fromPV)),weight=wgt)