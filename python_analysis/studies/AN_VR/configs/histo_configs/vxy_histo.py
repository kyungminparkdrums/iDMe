import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    h.make("sel_vtx_mindxy",'dxy')
    
    h.make("sel_vtx_vxy1",'vxy1')
    h.make("sel_vtx_vxy10",'vxy10')
    h.make("sel_vtx_vxy100",'vxy100')

    h.make("sel_vtx_vxy_fromPV1",'vxy_fromPV1')
    h.make("sel_vtx_vxy_fromPV10",'vxy_fromPV10')
    h.make("sel_vtx_vxy_fromPV100",'vxy_fromPV100')

    h.make("sel_vtx_vx1",'vx1')
    h.make("sel_vtx_vx10",'vx10')
    h.make("sel_vtx_vx100",'vx100')

    h.make("sel_vtx_vy1",'vy1')
    h.make("sel_vtx_vy10",'vy10')
    h.make("sel_vtx_vy100",'vy100')

    h.make("sel_vtx_vx_fromPV1",'vx_fromPV1')
    h.make("sel_vtx_vx_fromPV10",'vx_fromPV10')
    h.make("sel_vtx_vx_fromPV100",'vx_fromPV100')

    h.make("sel_vtx_vy_fromPV1",'vy_fromPV1')
    h.make("sel_vtx_vy_fromPV10",'vy_fromPV10')
    h.make("sel_vtx_vy_fromPV100",'vy_fromPV100')

    h.make("sel_vtx_CosThetaColl",'cosTheta')
    h.make("sel_vtx_CosThetaColl_fromPV",'cosTheta_fromPV')
    h.make("sel_vtx_CosThetaColl_fromPV_refit",'cosTheta_fromPV_refit')

    h.make("sel_vtx_vxy1_vs_CosThetaColl",'vxy1','cosTheta')
    h.make("sel_vtx_vxy10_vs_CosThetaColl",'vxy10','cosTheta')

    h.make("sel_vtx_vxy1_fromPV_vs_CosThetaColl_fromPV",'vxy_fromPV1','cosTheta_fromPV')
    h.make("sel_vtx_vxy10_fromPV_vs_CosThetaColl_fromPV",'vxy_fromPV10','cosTheta_fromPV')

    h.make("sel_vtx_vxy1_fromPV_vs_CosThetaColl_fromPV_refit",'vxy_fromPV1','cosTheta_fromPV_refit')
    h.make("sel_vtx_vxy10_fromPV_vs_CosThetaColl_fromPV_refit",'vxy_fromPV10','cosTheta_fromPV_refit')

    h.make("CosThetaColl_vs_CosThetaColl_fromPV",'cosTheta','cosTheta_fromPV')
    h.make("CosThetaColl_fromPV_vs_CosThetaColl_fromPV_refit",'cosTheta_fromPV','cosTheta_fromPV_refit')

    h.make("PVx",'pvx')
    h.make("PVy",'pvy')

    h.make("sel_vtx_sign",'vtx_sign')
    h.make("sel_vtx_sign_vs_sel_vtx_pt_ratio",'vtx_sign','ratio')
    h.make("sel_vtx_sign_vs_CosThetaColl_fromPV_refit",'vtx_sign','cosTheta_fromPV_refit')

    h.make("sel_vtx_pt_ratio",'ratio')
    h.make("sel_vtx_pt_ratio_vs_CosThetaColl_fromPV_refit",'ratio','cosTheta_fromPV_refit')
    
    h.make("lead_jet_phi_vs_lead_jet_eta","phi","eta")

    h.make("sel_vtx_dEta",'dEta')
    h.make("sel_vtx_dPhi",'dPhi')
    h.make("sel_vtx_log10dEtadPhi",'log10dEtadPhi')

    h.make("sel_vtx_log10dEtadPhi_vs_CosThetaColl_fromPV_refit", 'log10dEtadPhi', 'cosTheta_fromPV_refit')

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

    abs_dEta = np.abs(sel_vtx.e1.eta - sel_vtx.e2.eta)
    abs_dPhi = np.abs(sel_vtx.e1.phi - sel_vtx.e2.phi)
    
    ### FILLING HISTOGRAMS ###

    h.fill("sel_vtx_mindxy",dxy=sel_vtx.min_dxy,weight=wgt)
    
    h.fill("sel_vtx_vxy1",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy10",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy100",vxy=sel_vtx.vxy,weight=wgt)

    h.fill("sel_vtx_vxy_fromPV1",vxy_fromPV=vxy_fromPV,weight=wgt)
    h.fill("sel_vtx_vxy_fromPV10",vxy_fromPV=vxy_fromPV,weight=wgt)
    h.fill("sel_vtx_vxy_fromPV100",vxy_fromPV=vxy_fromPV,weight=wgt)

    h.fill("sel_vtx_vx1",vx=sel_vtx.vx,weight=wgt)
    h.fill("sel_vtx_vx10",vx=sel_vtx.vx,weight=wgt)
    h.fill("sel_vtx_vx100",vx=sel_vtx.vx,weight=wgt)

    h.fill("sel_vtx_vy1",vy=sel_vtx.vy,weight=wgt)
    h.fill("sel_vtx_vy10",vy=sel_vtx.vy,weight=wgt)
    h.fill("sel_vtx_vy100",vy=sel_vtx.vy,weight=wgt)

    h.fill("sel_vtx_vx_fromPV1",vx_fromPV=vx_fromPV,weight=wgt)
    h.fill("sel_vtx_vx_fromPV10",vx_fromPV=vx_fromPV,weight=wgt)
    h.fill("sel_vtx_vx_fromPV100",vx_fromPV=vx_fromPV,weight=wgt)

    h.fill("sel_vtx_vy_fromPV1",vy_fromPV=vy_fromPV,weight=wgt)
    h.fill("sel_vtx_vy_fromPV10",vy_fromPV=vy_fromPV,weight=wgt)
    h.fill("sel_vtx_vy_fromPV100",vy_fromPV=vy_fromPV,weight=wgt)

    h.fill("sel_vtx_CosThetaColl",cosTheta=sel_vtx.cos_collinear,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV",cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
    
    h.fill("sel_vtx_vxy1_vs_CosThetaColl",vxy=sel_vtx.vxy,cosTheta=sel_vtx.cos_collinear,weight=wgt)
    h.fill("sel_vtx_vxy10_vs_CosThetaColl",vxy=sel_vtx.vxy,cosTheta=sel_vtx.cos_collinear,weight=wgt)

    h.fill("sel_vtx_vxy1_fromPV_vs_CosThetaColl_fromPV",vxy_fromPV=vxy_fromPV,cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)
    h.fill("sel_vtx_vxy10_fromPV_vs_CosThetaColl_fromPV",vxy_fromPV=vxy_fromPV,cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)

    h.fill("sel_vtx_vxy1_fromPV_vs_CosThetaColl_fromPV_refit",\
           vxy_fromPV=vxy_fromPV,cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
    h.fill("sel_vtx_vxy10_fromPV_vs_CosThetaColl_fromPV_refit",\
           vxy_fromPV=vxy_fromPV,cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)

    h.fill("CosThetaColl_vs_CosThetaColl_fromPV",cosTheta=sel_vtx.cos_collinear,cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)        
    h.fill("CosThetaColl_fromPV_vs_CosThetaColl_fromPV_refit",\
            cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
    
    h.fill("PVx",pvx=events.PV.x,weight=wgt)
    h.fill("PVy",pvy=events.PV.y,weight=wgt)

    h.fill("sel_vtx_sign",sign=sel_vtx.sign,weight=wgt)
    h.fill("sel_vtx_sign_vs_sel_vtx_pt_ratio", sign=sel_vtx.sign, ratio=np.minimum(sel_vtx.e1.pt, sel_vtx.e2.pt)/np.maximum(sel_vtx.e1.pt, sel_vtx.e2.pt),weight=wgt)
    h.fill("sel_vtx_sign_vs_CosThetaColl_fromPV_refit", sign=sel_vtx.sign, cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)

    h.fill("sel_vtx_pt_ratio",ratio=np.minimum(sel_vtx.e1.pt, sel_vtx.e2.pt)/np.maximum(sel_vtx.e1.pt, sel_vtx.e2.pt),weight=wgt)
    h.fill("sel_vtx_pt_ratio_vs_CosThetaColl_fromPV_refit",\
           ratio=np.minimum(sel_vtx.e1.pt, sel_vtx.e2.pt)/np.maximum(sel_vtx.e1.pt, sel_vtx.e2.pt),\
           cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
    
    h.fill("lead_jet_phi_vs_lead_jet_eta",phi=events.PFJet.phi[:,0],eta=events.PFJet.eta[:,0],weight=wgt)

    h.fill("sel_vtx_dEta",dEta=abs_dEta,weight=wgt)
    h.fill("sel_vtx_dPhi",dPhi=abs_dPhi,weight=wgt)
    h.fill("sel_vtx_log10dEtadPhi",log10dEtadPhi=np.log10(abs_dEta/abs_dPhi),weight=wgt)

    h.fill("sel_vtx_log10dEtadPhi_vs_CosThetaColl_fromPV_refit", log10dEtadPhi=np.log10(abs_dEta/abs_dPhi), \
           cosTheta_fromPV_refit = sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
    
    #h.fill("bdtscore",bdtscore=events.BDTScore,weight=wgt)