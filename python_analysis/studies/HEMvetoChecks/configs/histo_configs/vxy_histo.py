import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    h.make("sel_vtx_vxy1",'vxy1')
    h.make("sel_vtx_vxy10",'vxy10')
    h.make("sel_vtx_vxy100",'vxy100')

    h.make("sel_vtx_vxyFromPV1",'vxyFromPV1')
    h.make("sel_vtx_vxyFromPV10",'vxyFromPV10')
    h.make("sel_vtx_vxyFromPV100",'vxyFromPV100')

    h.make("sel_vtx_vx1",'vx1')
    h.make("sel_vtx_vx10",'vx10')
    h.make("sel_vtx_vx100",'vx100')

    h.make("sel_vtx_vy1",'vy1')
    h.make("sel_vtx_vy10",'vy10')
    h.make("sel_vtx_vy100",'vy100')

    h.make("sel_vtx_vxFromPV1",'vxFromPV1')
    h.make("sel_vtx_vxFromPV10",'vxFromPV10')
    h.make("sel_vtx_vxFromPV100",'vxFromPV100')

    h.make("sel_vtx_vyFromPV1",'vyFromPV1')
    h.make("sel_vtx_vyFromPV10",'vyFromPV10')
    h.make("sel_vtx_vyFromPV100",'vyFromPV100')

    h.make("sel_vtx_CosThetaColl",'cosTheta')

    h.make("PVx",'pvx')
    h.make("PVy",'pvy')
    
    h.make("lead_jet_phi_vs_lead_jet_eta","phi","eta")

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

    vxFromPV = sel_vtx.vx - events.PV.x
    vyFromPV = sel_vtx.vy - events.PV.y

    vxyFromPV = np.sqrt(vxFromPV * vxFromPV + vyFromPV * vyFromPV)
    
    ### FILLING HISTOGRAMS ###

    h.fill("sel_vtx_vxy1",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy10",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy100",vxy=sel_vtx.vxy,weight=wgt)

    h.fill("sel_vtx_vxyFromPV1",vxyFromPV=vxyFromPV,weight=wgt)
    h.fill("sel_vtx_vxyFromPV10",vxyFromPV=vxyFromPV,weight=wgt)
    h.fill("sel_vtx_vxyFromPV100",vxyFromPV=vxyFromPV,weight=wgt)

    h.fill("sel_vtx_vx1",vx=sel_vtx.vx,weight=wgt)
    h.make("sel_vtx_vx10",vx=sel_vtx.vx,weight=wgt)
    h.make("sel_vtx_vx100",vx=sel_vtx.vx,weight=wgt)

    h.make("sel_vtx_vy1",vy=sel_vtx.vy,weight=wgt)
    h.make("sel_vtx_vy10",vy=sel_vtx.vy,weight=wgt)
    h.make("sel_vtx_vy100",vy=sel_vtx.vy,weight=wgt)

    h.fill("sel_vtx_vxFromPV1",vxFromPV=vxFromPV,weight=wgt)
    h.fill("sel_vtx_vxFromPV10",vxFromPV=vxFromPV,weight=wgt)
    h.fill("sel_vtx_vxFromPV100",vxFromPV=vxFromPV,weight=wgt)

    h.fill("sel_vtx_vyFromPV1",vyFromPV=vyFromPV,weight=wgt)
    h.fill("sel_vtx_vyFromPV10",vyFromPV=vyFromPV,weight=wgt)
    h.fill("sel_vtx_vyFromPV100",vyFromPV=vyFromPV,weight=wgt)

    h.fill("sel_vtx_CosThetaColl",cosTheta=sel_vtx.cos_collinear,weight=wgt)

    h.fill("PVx",pvx=events.PV.x,weight=wgt)
    h.fill("PVy",pvy=events.PV.y,weight=wgt)
    
    h.fill("lead_jet_phi_vs_lead_jet_eta",phi=events.PFJet.phi[:,0],eta=events.PFJet.eta[:,0],weight=wgt)
    
    #h.fill("bdtscore",bdtscore=events.BDTScore,weight=wgt)