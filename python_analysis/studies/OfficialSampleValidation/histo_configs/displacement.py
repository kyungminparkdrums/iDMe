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

    h.make("sel_e1_dxy_fromPV",'dxy')
    h.make("sel_e2_dxy_fromPV",'dxy')
    
    # selected vertex
    h.make("sel_vtx_dR",'dR')
    h.make("sel_vtx_mindxy",'dxy')
    h.make("sel_vtx_vxy1",'vxy1')
    h.make("sel_vtx_vxy10",'vxy10')
    h.make("sel_vtx_vxy20",'vxy20')
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
    #h.make("sel_vtx_CosThetaColl",'cosTheta')
    h.make("sel_vtx_LxyCosThetaColl",'LxyCosTheta')
    h.make("sel_vtx_LxyCosThetaCollZoom",'LxyCosThetaZoom')
    h.make("sel_vtx_LxyCosThetaCollZoomZoom",'LxyCosThetaZoomZoom')
    h.make("sel_vtx_eleDphi",'abs_dphi')
    h.make("sel_vtx_maxMiniRelIso",'iso')
    h.make("sel_vtx_maxMiniRelIsoCorr",'iso')

    ## Refitting of vtx track ##
    h.make('sel_vtx_mass_refit','vtx_mass_refit')
    h.make('sel_vtx_mass_low_refit','mass_low_refit')
    h.make('sel_vtx_pt_refit','vtx_pt_refit')
    h.make('sel_vtx_eta_refit','eta_refit')
    h.make('sel_vtx_phi_refit','phi_refit')
    h.make('sel_vtx_dR_refit','dR_refit')

    h.make('sel_vtx_mass_refit_vs_sel_vtx_mass','vtx_mass_refit','vtx_mass')
    h.make('sel_vtx_mass_low_refit_vs_sel_vtx_mass_low','mass_low_refit','mass_low')
    h.make('sel_vtx_pt_refit_vs_sel_vtx_pt','vtx_pt_refit','vtx_pt')
    h.make('sel_vtx_eta_refit_vs_sel_vtx_eta','eta_refit','eta')
    h.make('sel_vtx_phi_refit_vs_sel_vtx_phi','phi_refit','phi')
    h.make('sel_vtx_dR_refit_vs_sel_vtx_dR','dR_refit','dR')

    h.make("sel_vtx_mass_vs_vxy","mass_low","vxy10")
    h.make("sel_vtx_mass_refit_vs_vxy","mass_low_refit","vxy10")

    ##
    h.make("sel_e1_phi_vs_METphi","phi","metphi")
    h.make("sel_e1_phi_vs_sel_e1_eta","phi","eta")
    h.make("sel_e1_phi_vs_sel_e1_pt","phi","ele_pt")

    h.make("sel_e2_phi_vs_METphi","phi","metphi")
    h.make("sel_e2_phi_vs_sel_e2_eta","phi","eta")
    h.make("sel_e2_phi_vs_sel_e2_pt","phi","ele_pt")
    
    h.make("sel_vtx_phi_vs_sel_vtx_eta","phi","eta")
    h.make("sel_vtx_phi_vs_sel_vtx_pt","phi","vtx_pt")
    h.make("sel_vtx_phi_vs_sel_vtx_dR","phi","dR")

    h.make("sel_vtx_phi_vs_sel_vtx_pt800","phi","vtx_pt800")
    h.make("sel_vtx_phi_refit_zoom_vs_sel_vtx_phi_zoom","phi_refit_zoomHEM","phi_zoomHEM")
    h.make("sel_e1_phi_zoom_vs_sel_e1_eta","phi_zoomHEM","eta")
    h.make("sel_e2_phi_zoom_vs_sel_e2_eta","phi_zoomHEM","eta")

    ##
    
    h.make("sel_vtx_phi_vs_METphi","phi","metphi")
    h.make("sel_vtx_phi_vs_lead_jet_phi","phi","jet_phi")

    h.make("sel_vtx_phi_refit_vs_METphi","phi_refit","metphi")
    h.make("sel_vtx_phi_refit_vs_lead_jet_phi","phi_refit","jet_phi")
    
    h.make("METphi_vs_lead_jet_phi","metphi","jet_phi")

    h.make("sel_vtx_phi_refit_over_phi","ratio")

    h.make("sel_e1_log10dxydz","log10dEtadPhi")
    h.make("sel_e2_log10dxydz","log10dEtadPhi")
    h.make("min_sel_log10dxydz","log10dEtadPhi")
    
    ## ##
    
    h.make("sel_vtx_vx_vs_vy","vx","vy")
    h.make("PVx_vs_PVy","pvx","pvy")
    h.make("nPV",'nPV')
    
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
    h.make('lead_jet_phi','phi')
    h.make("jetMETratio",'jetMETratio')
    
    h.make("lead_jet_phi_vs_lead_jet_eta","phi","eta")

    h.make("PVx",'pvx')
    h.make("PVy",'pvy')

    h.make("sel_vtx_CosThetaColl",'cosTheta')

    h.make('sel_vtx_ThetaColl_fromPV_refit_rad','theta_rad')
    h.make('sel_vtx_ThetaColl_fromPV_refit','theta')
    h.make('sel_vtx_ThetaColl_fromPV','theta')

    h.make("sel_vtx_CosThetaColl_fromPV",'cosTheta_fromPV')
    h.make("sel_vtx_CosThetaColl_fromPV_refit",'cosTheta_fromPV_refit')

    h.make("sel_vtx_GenCosThetaColl",'cosTheta')
    h.make("sel_vtx_GenThetaColl",'gen_theta')
    h.make("sel_vtx_GenCosThetaColl_vs_sel_vtx_CosThetaColl_fromPV","cosTheta","cosTheta_fromPV")
    h.make("sel_vtx_GenCosThetaColl_vs_sel_vtx_CosThetaColl_fromPV_refit","cosTheta","cosTheta_fromPV_refit")

    h.make("sel_vtx_GenThetaColl_vs_sel_vtx_ThetaColl_fromPV","gen_theta","theta")
    h.make("sel_vtx_GenThetaColl_vs_sel_vtx_ThetaColl_fromPV_refit","gen_theta","theta")
    h.make("sel_vtx_GenThetaColl_rad_vs_sel_vtx_ThetaColl_fromPV_refit_rad","gen_theta_rad","theta_rad")

    h.make("sel_vtx_vxy1_fromPV_scalarSub",'vxy1')
    h.make("sel_vtx_vxy1_fromPV_vectorSub",'vxy1')

    h.make("sel_e1_dxy_refit", 'dxy')
    
    h.make("sel_e1_dxy_manual",'dxy')
    h.make("sel_e1_dxy_manual_fromPV",'dxy')

    h.make("bdtscore",'bdtscore')
    h.make("bdtscore_vs_sel_vtx_CosThetaColl_fromPV_refit","bdtscore","cosTheta_fromPV_refit")
    h.make('bdtscore_vs_sel_vtx_sign',"bdtscore",'vtx_sign')
    h.make('sel_vtx_CosThetaColl_fromPV_refit_vs_sel_vtx_sign',"cosTheta_fromPV_refit",'vtx_sign')

    h.make("bdtscore_vs_sel_vtx_CosThetaColl_fromPV","bdtscore","cosTheta_fromPV")
    h.make("bdtscore_vs_sel_vtx_thetaColl_fromPV_refit",'bdtscore','theta')
    h.make("bdtscore_vs_sel_vtx_thetaColl_fromPV",'bdtscore','theta')
    h.make('bdtscore_vs_sel_vtx_thetaColl_fromPV_refit_rad','bdtscore','theta_rad')

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
    #
    h.fill("sel_vtx_dR",dR=sel_vtx.dR,weight=wgt)
    #h.fill("sel_vtx_mindxy",dxy=sel_vtx.min_dxy,weight=wgt)
    h.fill("sel_vtx_vxy1",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy10",vxy=sel_vtx.vxy,weight=wgt)
    h.fill("sel_vtx_vxy20",vxy=sel_vtx.vxy,weight=wgt)

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
    #h.fill('sel_vtx_type',vtype=sel_vtx.typ,weight=wgt)
    #h.fill("sel_vtx_minEleDrJ",dR=np.minimum(sel_vtx.e1.mindRj,sel_vtx.e2.mindRj),weight=wgt)
    #h.fill("sel_vtx_minEleDPhiJ",abs_dphi=np.minimum(sel_vtx.e1.mindPhiJ,sel_vtx.e2.mindPhiJ),weight=wgt)
    h.fill("sel_vtx_mass_low",mass_low=sel_vtx.m,weight=wgt)
    h.fill("sel_vtx_mindxy_low",mindxy_low=sel_vtx.min_dxy,weight=wgt)
    #h.fill("sel_vtx_sign_etaProd",sign_etaProd=ak.values_astype(np.sign(sel_vtx.e1.eta*sel_vtx.e2.eta),int),weight=wgt)
    #h.fill("sel_vtx_CosThetaColl",cosTheta=sel_vtx.cos_collinear,weight=wgt)
    #h.fill("sel_vtx_LxyCosThetaColl",LxyCosTheta=sel_vtx.projectedLxy,weight=wgt)
    #h.fill("sel_vtx_LxyCosThetaCollZoom",LxyCosThetaZoom=sel_vtx.projectedLxy,weight=wgt)
    #h.fill("sel_vtx_LxyCosThetaCollZoomZoom",LxyCosThetaZoomZoom=sel_vtx.projectedLxy,weight=wgt)
    #h.fill("sel_vtx_eleDphi",abs_dphi=sel_vtx.eleDphi,weight=wgt)
    #h.fill("sel_vtx_maxMiniRelIso",iso=np.maximum(sel_vtx.e1.miniRelIso,sel_vtx.e2.miniRelIso),weight=wgt)
    #h.fill("sel_vtx_maxMiniRelIsoCorr",iso=np.maximum(sel_vtx.e1.miniRelIsoEleCorr,sel_vtx.e2.miniRelIsoEleCorr),weight=wgt)

    ## Refit
    
    h.fill('sel_vtx_mass_refit',mass_refit=sel_vtx.refit_m,weight=wgt)
    h.fill('sel_vtx_mass_low_refit',mass_low_refit=sel_vtx.refit_m,weight=wgt)
    h.fill('sel_vtx_pt_refit',pt_refit=sel_vtx.refit_pt,weight=wgt)
    h.fill('sel_vtx_eta_refit',eta_refit=sel_vtx.refit_eta,weight=wgt)
    h.fill('sel_vtx_phi_refit',phi_refit=sel_vtx.refit_phi,weight=wgt)
    h.fill('sel_vtx_dR_refit',dR_refit=sel_vtx.refit_dR,weight=wgt)