import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    # electron
    h.make("sel_e1_pt",'ele_pt')
    h.make("sel_e1_eta",'eta')
    h.make("sel_e1_phi",'phi')
    h.make("sel_e1_dxy",'dxy')
    h.make("sel_e1_dxySignif",'dxy_signif')
    h.make("sel_e1_vxy1",'vxy1')
    h.make("sel_e1_vz",'vz')

    h.make("sel_e2_pt",'ele_pt')
    h.make("sel_e2_eta",'eta')
    h.make("sel_e2_phi",'phi')
    h.make("sel_e2_dxy",'dxy')
    h.make("sel_e2_dxySignif",'dxy_signif')
    h.make("sel_e2_vxy1",'vxy1')
    h.make("sel_e2_vz",'vz')
    
    # selected vertex
    h.make("sel_vtx_dR",'dR')
    h.make("sel_vtx_mindxy",'dxy')
    h.make("sel_vtx_vxy1",'vxy1')
    h.make("sel_vtx_vxy10",'vxy10')
    h.make("sel_vtx_vxy100",'vxy100')
    h.make("sel_vtx_leadpT",'ele_pt')
    h.make("sel_vtx_METdPhi",'abs_dphi')
    h.make("sel_vtx_mindRj",'dR')
    h.make("sel_vtx_chi2",'vtx_chi2')
    h.make('sel_vtx_mass','vtx_mass')
    h.make('sel_vtx_mindPhiJ','abs_dphi')
    h.make('sel_vtx_sign','vtx_sign')
    h.make('sel_vtx_pt','vtx_pt')
    h.make('sel_vtx_eta','eta')
    h.make('sel_vtx_phi','phi')
    h.make("sel_vtx_type",'vtx_type')
    h.make("sel_vtx_minEleDrJ",'dR')
    h.make("sel_vtx_minEleDPhiJ",'abs_dphi')
    h.make("sel_vtx_mass_low",'mass_low')
    h.make("sel_vtx_mindxy_low",'mindxy_low')
    h.make("sel_vtx_sign_etaProd",'sign_etaProd')
    h.make("sel_vtx_CosThetaColl",'cosTheta')
    h.make("sel_vtx_LxyCosThetaColl",'LxyCosTheta')
    h.make("sel_vtx_LxyCosThetaCollZoom",'LxyCosThetaZoom')
    h.make("sel_vtx_LxyCosThetaCollZoomZoom",'LxyCosThetaZoomZoom')
    h.make("sel_vtx_eleDphi",'abs_dphi')
    h.make("sel_vtx_maxMiniRelIso",'iso')
    h.make("sel_vtx_maxMiniRelIsoCorr",'iso')
    
    h.make("sel_vtx_vx_vs_vy","vx","vy")
    h.make("sel_vtx_phi_vs_METphi","phi","metphi")

    h.make("PVx_vs_PVy","pvx","pvy")
    
    #h.make("sel_vtx_mass_vs_vxy","mass_low","vxy10")
    #h.make("sel_vtx_mass_vs_vxy20","mass_low","vxy20")
    #h.make("sel_vtx_vxyRes","vxyRes")
    #h.make("sel_vtx_vxySignif","vxySignif")
    #h.make("sel_vtx_vxy20_vs_vxyRes","vxy20","vxyRes")
    #h.make("sel_vtx_vxy20_vs_vxySignif","vxy20","vxySignif")
    
    # misc
    h.make("PFMET",'met')
    h.make("PFMET1000",'met1000')
    h.make("PFMETphi",'metphi')
    h.make("jetMETdPhi",'abs_dphi')
    h.make("minJetMETdPhi",'abs_dphi')
    h.make("nJets",'nJets')
    h.make('lead_jet_pt','jet_pt')
    h.make('lead_jet_eta','eta')
    h.make("jetMETratio",'jetMETratio')

    h.make("PVx",'pvx')
    h.make("PVy",'pvy')

    #h.make("bdtscore",'bdtscore')

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
    h.fill("sel_e1_eta",eta=sel_vtx.e1.eta,weight=wgt)
    h.fill("sel_e1_phi",phi=sel_vtx.e1.phi,weight=wgt)
    h.fill("sel_e1_dxy",dxy=sel_vtx.e1.dxy,weight=wgt)
    h.fill("sel_e1_vxy1",vxy=sel_vtx.e1.vxy,weight=wgt)
    h.fill("sel_e1_vz",vz=sel_vtx.e1.vz,weight=wgt)

    h.fill("sel_e2_pt",pt=sel_vtx.e2.pt,weight=wgt)
    h.fill("sel_e2_eta",eta=sel_vtx.e2.eta,weight=wgt)
    h.fill("sel_e2_phi",phi=sel_vtx.e2.phi,weight=wgt)
    h.fill("sel_e2_dxy",dxy=sel_vtx.e2.dxy,weight=wgt)
    h.fill("sel_e2_vxy1",vxy=sel_vtx.e2.vxy,weight=wgt)
    h.fill("sel_e2_vz",vz=sel_vtx.e2.vz,weight=wgt)
    
    #
    h.fill("sel_vtx_dR",dR=sel_vtx.dR,weight=wgt)
    h.fill("sel_vtx_mindxy",dxy=sel_vtx.min_dxy,weight=wgt)
    h.fill("sel_vtx_vxy1",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy10",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy100",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_leadpT",pt=np.maximum(sel_vtx.e1.pt,sel_vtx.e2.pt),weight=wgt)
    h.fill("sel_vtx_METdPhi",abs_dphi=np.abs(sel_vtx.METdPhi),weight=wgt)
    h.fill('sel_vtx_mindRj',dR=sel_vtx.mindRj,weight=wgt)
    h.fill('sel_vtx_chi2',chi2=sel_vtx.reduced_chi2,weight=wgt)
    h.fill('sel_vtx_mass',mass=sel_vtx.m,weight=wgt)
    h.fill('sel_vtx_mindPhiJ',abs_dphi=np.abs(sel_vtx.mindPhiJ),weight=wgt)
    h.fill('sel_vtx_sign',sign=sel_vtx.sign,weight=wgt)
    h.fill('sel_vtx_pt',pt=sel_vtx.pt,weight=wgt)
    h.fill('sel_vtx_eta',eta=sel_vtx.eta,weight=wgt)
    h.fill('sel_vtx_phi',phi=sel_vtx.phi,weight=wgt)
    h.fill('sel_vtx_type',vtype=sel_vtx.typ,weight=wgt)
    h.fill("sel_vtx_minEleDrJ",dR=np.minimum(sel_vtx.e1.mindRj,sel_vtx.e2.mindRj),weight=wgt)
    h.fill("sel_vtx_minEleDPhiJ",abs_dphi=np.minimum(sel_vtx.e1.mindPhiJ,sel_vtx.e2.mindPhiJ),weight=wgt)
    h.fill("sel_vtx_mass_low",mass_low=sel_vtx.m,weight=wgt)
    h.fill("sel_vtx_mindxy_low",mindxy_low=sel_vtx.min_dxy,weight=wgt)
    h.fill("sel_vtx_sign_etaProd",sign_etaProd=ak.values_astype(np.sign(sel_vtx.e1.eta*sel_vtx.e2.eta),int),weight=wgt)
    h.fill("sel_vtx_CosThetaColl",cosTheta=sel_vtx.cos_collinear,weight=wgt)
    h.fill("sel_vtx_LxyCosThetaColl",LxyCosTheta=sel_vtx.projectedLxy,weight=wgt)
    h.fill("sel_vtx_LxyCosThetaCollZoom",LxyCosThetaZoom=sel_vtx.projectedLxy,weight=wgt)
    h.fill("sel_vtx_LxyCosThetaCollZoomZoom",LxyCosThetaZoomZoom=sel_vtx.projectedLxy,weight=wgt)
    h.fill("sel_vtx_eleDphi",abs_dphi=sel_vtx.eleDphi,weight=wgt)
    h.fill("sel_vtx_maxMiniRelIso",iso=np.maximum(sel_vtx.e1.miniRelIso,sel_vtx.e2.miniRelIso),weight=wgt)
    h.fill("sel_vtx_maxMiniRelIsoCorr",iso=np.maximum(sel_vtx.e1.miniRelIsoEleCorr,sel_vtx.e2.miniRelIsoEleCorr),weight=wgt)
    
    h.fill("sel_vtx_vx_vs_vy",vx=sel_vtx.vx,vy=sel_vtx.vy,weight=wgt)
    h.fill("sel_vtx_phi_vs_METphi",phi=sel_vtx.phi,metphi=events.PFMET.phi,weight=wgt)
    h.fill("PVx_vs_PVy",pvx=events.PV.x,pvy=events.PV.y,weight=wgt)
    
    #h.fill("sel_vtx_mass_vs_vxy",mass_low=sel_vtx.m,vxy=sel_vtx.vxy,weight=wgt)
    #h.fill("sel_vtx_mass_vs_vxy20",mass_low=sel_vtx.m,vxy=sel_vtx.vxy,weight=wgt)
    #h.fill("sel_vtx_vxyRes",vxyRes=sel_vtx.sigmavxy,weight=wgt)
    #h.fill("sel_vtx_vxySignif",vxySignif=sel_vtx.vxy/sel_vtx.sigmavxy,weight=wgt)
    #h.fill("sel_vtx_vxy20_vs_vxyRes",vxy=sel_vtx.vxy,vxyRes=sel_vtx.sigmavxy,weight=wgt)
    #h.fill("sel_vtx_vxy20_vs_vxySignif",vxy=sel_vtx.vxy,vxySignif=sel_vtx.vxy/sel_vtx.sigmavxy,weight=wgt)

    #
    h.fill("PFMET",met=events.PFMET.pt,weight=wgt)
    h.fill("PFMET1000",met=events.PFMET.pt,weight=wgt)
    h.fill("PFMETphi",metphi=events.PFMET.phi,weight=wgt)
    h.fill("jetMETdPhi",abs_dphi=np.abs(events.PFJet.METdPhi[:,0]),weight=wgt)
    h.fill("minJetMETdPhi",abs_dphi=ak.min(np.abs(events.PFJet.METdPhi),axis=1),weight=wgt)
    h.fill("nJets",nJets=ak.count(events.PFJet.pt,axis=1),weight=wgt)
    h.fill("lead_jet_pt",jet_pt=events.PFJet.pt[:,0],weight=wgt)
    h.fill("lead_jet_eta",eta=events.PFJet.eta[:,0],weight=wgt)
    h.fill("jetMETratio",jetMETratio=events.PFJet.pt[:,0]/events.PFMET.pt,weight=wgt)

    h.fill("PVx",pvx=events.PV.x,weight=wgt)
    h.fill("PVy",pvy=events.PV.y,weight=wgt)

    #h.fill("bdtscore",bdtscore=events.BDTScore,weight=wgt)