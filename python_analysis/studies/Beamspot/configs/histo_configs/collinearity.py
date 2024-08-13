import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    h.make("sel_vtx_CosThetaColl",'cosTheta')
    h.make("sel_vtx_CosThetaColl_fromPV",'cosTheta_fromPV')
    #h.make("sel_vtx_CosThetaColl_fromPV_refit",'cosTheta_fromPV_refit')

    h.make("sel_vtx_vxy1_vs_CosThetaColl",'vxy1','cosTheta')
    h.make("sel_vtx_vxy10_vs_CosThetaColl",'vxy10','cosTheta')

    h.make("sel_vtx_vxy1_fromPV_vs_CosThetaColl_fromPV",'vxy_fromPV1','cosTheta_fromPV')
    h.make("sel_vtx_vxy10_fromPV_vs_CosThetaColl_fromPV",'vxy_fromPV10','cosTheta_fromPV')

    #h.make("sel_vtx_vxy1_fromPV_vs_CosThetaColl_fromPV_refit",'vxy_fromPV1','cosTheta_fromPV_refit')
    #h.make("sel_vtx_vxy10_fromPV_vs_CosThetaColl_fromPV_refit",'vxy_fromPV10','cosTheta_fromPV_refit')

    h.make("CosThetaColl_vs_CosThetaColl_fromPV",'cosTheta','cosTheta_fromPV')
    #h.make("CosThetaColl_fromPV_vs_CosThetaColl_fromPV_refit",'cosTheta_fromPV','cosTheta_fromPV_refit')

    h.make("PVx",'pvx')
    h.make("PVy",'pvy')

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

    vx_fromPV = sel_vtx.vx - events.PV.x
    vy_fromPV = sel_vtx.vy - events.PV.y

    vxy_fromPV = np.sqrt(vx_fromPV * vx_fromPV + vy_fromPV * vy_fromPV)
    
    ### FILLING HISTOGRAMS ###

    h.fill("sel_vtx_CosThetaColl",cosTheta=sel_vtx.cos_collinear,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV",cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)
    #h.fill("sel_vtx_CosThetaColl_fromPV_refit",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
    
    h.fill("sel_vtx_vxy1_vs_CosThetaColl",vxy=sel_vtx.vxy,cosTheta=sel_vtx.cos_collinear,weight=wgt)
    h.fill("sel_vtx_vxy10_vs_CosThetaColl",vxy=sel_vtx.vxy,cosTheta=sel_vtx.cos_collinear,weight=wgt)

    h.fill("sel_vtx_vxy1_fromPV_vs_CosThetaColl_fromPV",vxy_fromPV=vxy_fromPV,cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)
    h.fill("sel_vtx_vxy10_fromPV_vs_CosThetaColl_fromPV",vxy_fromPV=vxy_fromPV,cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)

    #h.fill("sel_vtx_vxy1_fromPV_vs_CosThetaColl_fromPV_refit",\
    #       vxy_fromPV=vxy_fromPV,cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
    #h.fill("sel_vtx_vxy10_fromPV_vs_CosThetaColl_fromPV_refit",\
    #       vxy_fromPV=vxy_fromPV,cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)

    h.fill("CosThetaColl_vs_CosThetaColl_fromPV",cosTheta=sel_vtx.cos_collinear,cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)        
    #h.fill("CosThetaColl_fromPV_vs_CosThetaColl_fromPV_refit",\
    #        cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
    
    h.fill("PVx",pvx=events.PV.x,weight=wgt)
    h.fill("PVy",pvy=events.PV.y,weight=wgt)