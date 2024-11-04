import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    h.make("sel_e1_dxy", "dxy")
    h.make("sel_e1_trkVtxVx", "trkVtxVx")
    h.make("sel_e1_trkVtxVy", "trkVtxVy")
    h.make("sel_e1_trkVtxPx", "trkVtxPx")
    h.make("sel_e1_trkVtxPy", "trkVtxPy")
    h.make("sel_e1_trkVtxPt", "trkVtxPt")
    h.make("sel_e1_trkVtxPVposx", "trkVtxPVposx")
    h.make("sel_e1_trkVtxPVposy", "trkVtxPVposy")
    h.make("sel_e1_trkVtxDxyFromPV", "dxy")
    h.make("sel_e1_dxy_refit", 'dxy')

    h.make("sel_e2_dxy", "dxy")
    h.make("sel_e2_trkVtxVx", "trkVtxVx")
    h.make("sel_e2_trkVtxVy", "trkVtxVy")
    h.make("sel_e2_trkVtxPx", "trkVtxPx")
    h.make("sel_e2_trkVtxPy", "trkVtxPy")
    h.make("sel_e2_trkVtxPt", "trkVtxPt")
    h.make("sel_e2_trkVtxPVposx", "trkVtxPVposx")
    h.make("sel_e2_trkVtxPVposy", "trkVtxPVposy")
    h.make("sel_e2_trkVtxDxyFromPV", "dxy")
    h.make("sel_e2_dxy_refit", 'dxy')

    h.make("sel_e_lead_dxy", "dxy")
    h.make("sel_e_sub_dxy", "dxy")
    h.make("sel_e_lead_dxy_refit", "dxy")
    h.make("sel_e_sub_dxy_refit", "dxy")
    
    h.make("PVx",'pvx')
    h.make("PVy",'pvy')

    h.make("sel_vtx_vxy1",'vxy1')
    h.make("sel_vtx_vxy10",'vxy10')
    h.make("sel_vtx_vxy100",'vxy100')

    h.make("sel_vtx_vxy1_fromPV_vectorSub",'vxy1')

    h.make("sel_vtx_CosThetaColl",'cosTheta')
    h.make("sel_vtx_CosThetaColl_fromPV",'cosTheta_fromPV')
    h.make("sel_vtx_CosThetaColl_fromPV_refit",'cosTheta_fromPV_refit')

    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e_lead_dxy",'cosTheta_fromPV_refit','dxy_zoom')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e_sub_dxy",'cosTheta_fromPV_refit','dxy_zoom')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e_lead_dxy_refit",'cosTheta_fromPV_refit','dxy_zoom')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e_sub_dxy_refit",'cosTheta_fromPV_refit','dxy_zoom')

    h.make("sel_vtx_min_dxy", "dxy")
    h.make("sel_vtx_min_dxy_refit", "dxy")

    h.make("sel_e1_trkVtxVxMinusBSx", "trkVtxVx")
    h.make("sel_e1_trkVtxVyMinusBSy", "trkVtxVy")

    h.make("sel_e1_trkVtxVxMinusPVx", "trkVtxVx")
    h.make("sel_e1_trkVtxVyMinusPVy", "trkVtxVy")

    h.make("PFMET",'met')
    h.make("PFMET1000",'met1000')
    h.make("PFMETphi",'metphi')
    h.make("jetMETdPhi",'abs_dphi')
    h.make("minJetMETdPhi",'abs_dphi')
    h.make("nJets",'nJets')
    h.make('lead_jet_pt','jet_pt')
    h.make('lead_jet_eta','eta')
    h.make('lead_jet_phi','phi')
    h.make("jetMETratio",'jetMETratio')

    h.make("lead_jet_phi_vs_lead_jet_eta","phi","eta")
    h.make("sel_vtx_phi_vs_METphi","phi","metphi")
    h.make("METphi_vs_lead_jet_phi","metphi","jet_phi")

    h.make('sel_vtx_mass_refit','vtx_mass_refit')
    h.make('sel_vtx_dR_refit','dR_refit')
    
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
    h.fill("sel_e1_dxy",dxy=sel_vtx.e1.dxy,weight=wgt)
    h.fill("sel_e1_trkVtxVx",trkVtxVx=sel_vtx.e1.trkVtxVx,weight=wgt)
    h.fill("sel_e1_trkVtxVy",trkVtxVy=sel_vtx.e1.trkVtxVy,weight=wgt)
    h.fill("sel_e1_trkVtxPx",trkVtxPx=sel_vtx.e1.trkVtxPx,weight=wgt)
    h.fill("sel_e1_trkVtxPy",trkVtxPy=sel_vtx.e1.trkVtxPy,weight=wgt)
    h.fill("sel_e1_trkVtxPt",trkVtxPt=sel_vtx.e1.trkVtxPt,weight=wgt)
    h.fill("sel_e1_trkVtxPVposx",trkVtxPVposx=sel_vtx.e1.trkPVposx,weight=wgt)
    h.fill("sel_e1_trkVtxPVposy",trkVtxPVposy=sel_vtx.e1.trkPVposy,weight=wgt)
    h.fill("sel_e1_trkVtxDxyFromPV",dxy=sel_vtx.e1.trkDxyFromPV,weight=wgt)
    h.fill("sel_e1_dxy_refit", dxy=np.abs(sel_vtx.e1.refit_dxy), weight=wgt)

    h.fill("sel_e2_dxy",dxy=sel_vtx.e2.dxy,weight=wgt)
    h.fill("sel_e2_trkVtxVx",trkVtxVx=sel_vtx.e2.trkVtxVx,weight=wgt)
    h.fill("sel_e2_trkVtxVy",trkVtxVy=sel_vtx.e2.trkVtxVy,weight=wgt)
    h.fill("sel_e2_trkVtxPx",trkVtxPx=sel_vtx.e2.trkVtxPx,weight=wgt)
    h.fill("sel_e2_trkVtxPy",trkVtxPy=sel_vtx.e2.trkVtxPy,weight=wgt)
    h.fill("sel_e2_trkVtxPt",trkVtxPt=sel_vtx.e2.trkVtxPt,weight=wgt)
    h.fill("sel_e2_trkVtxPVposx",trkVtxPVposx=sel_vtx.e2.trkPVposx,weight=wgt)
    h.fill("sel_e2_trkVtxPVposy",trkVtxPVposy=sel_vtx.e2.trkPVposy,weight=wgt)
    h.fill("sel_e2_trkVtxDxyFromPV",dxy=sel_vtx.e2.trkDxyFromPV,weight=wgt)
    h.fill("sel_e2_dxy_refit", dxy=np.abs(sel_vtx.e2.refit_dxy), weight=wgt)

    h.fill("sel_e_lead_dxy", dxy=np.maximum(sel_vtx.e1.dxy, sel_vtx.e2.dxy), weight=wgt)
    h.fill("sel_e_sub_dxy", dxy=np.minimum(sel_vtx.e1.dxy, sel_vtx.e2.dxy), weight=wgt)
    h.fill("sel_e_lead_dxy_refit", dxy=np.maximum(np.abs(sel_vtx.e1.refit_dxy), np.abs(sel_vtx.e2.refit_dxy)), weight=wgt)
    h.fill("sel_e_sub_dxy_refit", dxy=np.minimum(np.abs(sel_vtx.e1.refit_dxy), np.abs(sel_vtx.e2.refit_dxy)), weight=wgt)
    
    h.fill("PVx",pvx=events.PV.x,weight=wgt)
    h.fill("PVy",pvy=events.PV.y,weight=wgt)

    h.fill("sel_vtx_vxy1",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy10",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy100",vxy=sel_vtx.vxy,weight=wgt)

    h.fill("sel_vtx_vxy1_fromPV_vectorSub", vxy=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2), weight=wgt)

    h.fill("sel_vtx_CosThetaColl",cosTheta=sel_vtx.cos_collinear,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV",cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)

    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e_lead_dxy",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dxy=np.maximum(sel_vtx.e1.dxy, sel_vtx.e2.dxy), weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e_sub_dxy",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dxy=np.minimum(sel_vtx.e1.dxy, sel_vtx.e2.dxy), weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e_lead_dxy_refit",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dxy=np.maximum(np.abs(sel_vtx.e1.refit_dxy), np.abs(sel_vtx.e2.refit_dxy)), weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e_sub_dxy_refit",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dxy=np.minimum(np.abs(sel_vtx.e1.refit_dxy), np.abs(sel_vtx.e2.refit_dxy)), weight=wgt)

    h.fill("sel_vtx_min_dxy", dxy=np.minimum(sel_vtx.e1.dxy, sel_vtx.e2.dxy), weight=wgt)
    h.fill("sel_vtx_min_dxy_refit", dxy=np.minimum(np.abs(sel_vtx.e1.refit_dxy), np.abs(sel_vtx.e2.refit_dxy)), weight=wgt)

    h.fill("sel_e1_trkVtxVxMinusBSx", trkVtxVx=sel_vtx.e1.trkVtxVx-sel_vtx.e1.trkPVposx, weight=wgt)
    h.fill("sel_e1_trkVtxVyMinusBSy", trkVtxVy=sel_vtx.e1.trkVtxVy-sel_vtx.e1.trkPVposy, weight=wgt)

    h.fill("sel_e1_trkVtxVxMinusPVx", trkVtxVx=sel_vtx.e1.trkVtxVx-events.PV.x, weight=wgt)
    h.fill("sel_e1_trkVtxVyMinusPVy", trkVtxVy=sel_vtx.e1.trkVtxVy-events.PV.y, weight=wgt)

    h.fill("PFMET",met=events.PFMET.pt,weight=wgt)
    h.fill("PFMET1000",met=events.PFMET.pt,weight=wgt)
    h.fill("PFMETphi",metphi=events.PFMET.phi,weight=wgt)
    h.fill("jetMETdPhi",abs_dphi=np.abs(events.PFJet.METdPhi[:,0]),weight=wgt)
    h.fill("minJetMETdPhi",abs_dphi=ak.min(np.abs(events.PFJet.METdPhi),axis=1),weight=wgt)
    h.fill("nJets",nJets=ak.count(events.PFJet.pt,axis=1),weight=wgt)
    h.fill("lead_jet_pt",jet_pt=events.PFJet.pt[:,0],weight=wgt)
    h.fill("lead_jet_eta",eta=events.PFJet.eta[:,0],weight=wgt)
    h.fill("lead_jet_phi",phi=events.PFJet.phi[:,0],weight=wgt)
    h.fill("jetMETratio",jetMETratio=events.PFJet.pt[:,0]/events.PFMET.pt,weight=wgt)

    h.fill("lead_jet_phi_vs_lead_jet_eta",phi=events.PFJet.phi[:,0],eta=events.PFJet.eta[:,0],weight=wgt)
    h.fill("sel_vtx_phi_vs_METphi",phi=sel_vtx.phi,metphi=events.PFMET.phi,weight=wgt)
    h.fill("METphi_vs_lead_jet_phi",metphi=events.PFMET.phi,jet_phi=events.PFJet.phi[:,0],weight=wgt)

    h.fill('sel_vtx_mass_refit',mass_refit=sel_vtx.refit_m,weight=wgt)
    h.fill('sel_vtx_dR_refit',dR_refit=sel_vtx.refit_dR,weight=wgt)
    