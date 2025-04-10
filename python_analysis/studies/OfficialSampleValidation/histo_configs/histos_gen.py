import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    h.make("genEle_pt","ele_pt")
    h.make("genEle_eta","eta")
    h.make("genEle_phi","phi")
    h.make("genEle_vxy","vxy20")

    h.make("genPos_pt","ele_pt")
    h.make("genPos_eta","eta")
    h.make("genPos_phi","phi")
    h.make("genPos_vxy","vxy20")
    
    h.make("genEE_pt","vtx_pt")
    h.make("genEE_eta","eta")
    h.make("genEE_phi","phi")
    h.make("genEE_mass","vtx_mass")
    h.make("genEE_dR","dR")
    h.make("genEE_METdPhi","abs_dphi")
    h.make("genEE_vxy","vxy20")

    h.make("gen_JetMETdPhi_all","abs_dphi")
    h.make("gen_JetMETdPhi_lead","abs_dphi")
    h.make("gen_jetMETratio_lead","jetMETratio")
    h.make("genLeadMETPt","met1000")
    h.make("genLeadMETPhi","phi")

    h.make("genJetPt_all","jet_pt")
    h.make("genJetPt_lead","jet_pt")
    h.make("genJetEta_all","eta")
    h.make("genJetEta_lead","eta")
    h.make("genJetPhi_all","phi")
    h.make("genJetPhi_lead","phi")

    h.make("nGenJet_all",'nJets')
    h.make("nGenJet_30","nJets")
    
    h.make("sel_vtx_GenCosThetaColl","cosTheta")
    h.make("sel_vtx_GenThetaColl",'cosTheta')

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

    if info["type"] == "signal":
        h.fill("genEle_pt",pt=events.GenEle.pt,weight=wgt)
        h.fill("genEle_eta",eta=events.GenEle.eta,weight=wgt)
        h.fill("genEle_phi",phi=events.GenEle.phi,weight=wgt)
        h.fill("genEle_vxy",vxy=events.GenEle.vxy,weight=wgt)

        h.fill("genPos_pt",pt=events.GenPos.pt,weight=wgt)
        h.fill("genPos_eta",eta=events.GenPos.eta,weight=wgt)
        h.fill("genPos_phi",phi=events.GenPos.phi,weight=wgt)
        h.fill("genPos_vxy",vxy=events.GenPos.vxy,weight=wgt)

        h.fill("genEE_pt",pt=events.genEE.pt,weight=wgt)
        h.fill("genEE_eta",eta=events.genEE.eta,weight=wgt)
        h.fill("genEE_phi",phi=events.genEE.phi,weight=wgt)
        h.fill("genEE_mass",mass=events.genEE.mass,weight=wgt)
        h.fill("genEE_dR",dR=events.genEE.dr,weight=wgt)
        h.fill("genEE_METdPhi",abs_dphi=events.genEE.METdPhi,weight=wgt)
        h.fill("genEE_vxy",vxy=events.genEE.vxy,weight=wgt)

        #h.make("gen_JetMETdPhi_all","abs_dphi")
        h.fill("gen_JetMETdPhi_lead",abs_dphi=np.abs(events.GenJet.phi[:,0]-events.GenMET.phi),weight=wgt)
        h.fill("gen_jetMETratio_lead",jetMETratio=np.abs(events.GenMET.pt/events.GenJet.pt[:,0]),weight=wgt)
        h.fill("genLeadMETPt",met=events.GenMET.pt,weight=wgt)
        h.fill("genLeadMETPhi",phi=events.GenMET.phi,weight=wgt)
    
        #h.make("genJetPt_all","jet_pt")
        h.fill("genJetPt_lead",jet_pt=events.GenJet.pt[:,0],weight=wgt)
        #h.make("genJetEta_all","eta")
        h.fill("genJetEta_lead",eta=events.GenJet.eta[:,0],weight=wgt)
        #h.make("genJetPhi_all","phi")
        h.fill("genJetPhi_lead",phi=events.GenJet.phi[:,0],weight=wgt)
        
        h.fill("nGenJet_all",nJets=ak.num(events.GenJet.pt, axis=1),weight=wgt)
        h.fill("nGenJet_30",nJets=ak.count_nonzero(events.GenJet.pt[events.GenJet.pt>30], axis=1),weight=wgt)
        
        h.fill("sel_vtx_GenCosThetaColl",cosTheta=sel_vtx.gen_cos_collinear_fromPV,weight=wgt)
        h.fill("sel_vtx_GenThetaColl",cosTheta=np.degrees(np.arccos(sel_vtx.gen_cos_collinear_fromPV)),weight=wgt)
        