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
    h.make("sel_e1_dxy_fromPV",'dxy')
    
    h.make("sel_e1_etaErr",'eta_err')
    h.make("sel_e2_etaErr",'eta_err')
    h.make("sel_e1_phiErr",'phi_err')
    h.make("sel_e2_phiErr",'phi_err')
    
    # electron track stuff
    h.make("sel_e1_dz",'dz')
    h.make("sel_e1_dzErr",'dz_err')
    h.make("sel_e1_dzSignif",'dz_signif')
    h.make("sel_e1_chi2",'trk_chi2')
    h.make("sel_e1_prob",'prob')
    h.make("sel_e1_angularRes",'angular_res')
    h.make("sel_e1_trkiso",'iso')
    h.make("sel_e1_trk_reliso",'iso')

    h.make("sel_e2_dz",'dz')
    h.make("sel_e2_dzErr",'dz_err')
    h.make("sel_e2_dzSignif",'dz_signif')
    h.make("sel_e2_chi2",'trk_chi2')
    h.make("sel_e2_prob",'prob')
    h.make("sel_e2_angularRes",'angular_res')
    h.make("sel_e2_trkiso",'iso')
    h.make("sel_e2_trk_reliso",'iso')
    
    h.make("sel_vtx_CosThetaColl",'cosTheta')

    h.make('sel_vtx_ThetaColl_fromPV_refit','theta')
    h.make('sel_vtx_ThetaColl_fromPV','theta')

    h.make("sel_vtx_CosThetaColl_fromPV",'cosTheta_fromPV')
    h.make("sel_vtx_CosThetaColl_fromPV_refit",'cosTheta_fromPV_refit')

    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_etaErr",'cosTheta_fromPV_refit','eta_err')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_etaErr",'cosTheta_fromPV_refit','eta_err')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_phiErr",'cosTheta_fromPV_refit','phi_err')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_phiErr",'cosTheta_fromPV_refit','phi_err')

    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_dz",'cosTheta_fromPV_refit','dz')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_dzErr",'cosTheta_fromPV_refit','dz_err')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_dzSignif",'cosTheta_fromPV_refit','dz_signif')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_chi2",'cosTheta_fromPV_refit','trk_chi2')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_prob",'cosTheta_fromPV_refit','prob')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_angularRes",'cosTheta_fromPV_refit','angular_res')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_trkiso",'cosTheta_fromPV_refit','iso')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_trk_reliso",'cosTheta_fromPV_refit','iso')

    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_dz",'cosTheta_fromPV_refit','dz')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_dzErr",'cosTheta_fromPV_refit','dz_err')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_dzSignif",'cosTheta_fromPV_refit','dz_signif')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_chi2",'cosTheta_fromPV_refit','trk_chi2')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_prob",'cosTheta_fromPV_refit','prob')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_angularRes",'cosTheta_fromPV_refit','angular_res')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_trkiso",'cosTheta_fromPV_refit','iso')
    h.make("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_trk_reliso",'cosTheta_fromPV_refit','iso')


    h.make('sel_vtx_mass_refit','vtx_mass_refit')
    h.make('sel_vtx_dR_refit','dR_refit')
    
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

    h.fill("sel_e1_dxySignif",signif=sel_vtx.e1.dxy/sel_vtx.e1.dxyErr,weight=wgt)
    h.fill("sel_e2_dxySignif",signif=sel_vtx.e2.dxy/sel_vtx.e2.dxyErr,weight=wgt)
    
    h.fill("sel_e1_etaErr",eta_err=sel_vtx.e1.etaErr,weight=wgt)
    h.fill("sel_e2_etaErr",eta_err=sel_vtx.e2.etaErr,weight=wgt)
    h.fill("sel_e1_phiErr",phi_err=sel_vtx.e1.phiErr,weight=wgt)
    h.fill("sel_e2_phiErr",phi_err=sel_vtx.e2.phiErr,weight=wgt)
    
    h.fill("sel_e1_dz",dz=sel_vtx.e1.dz,weight=wgt)
    h.fill("sel_e1_dzErr",dz_err=sel_vtx.e1.dzErr,weight=wgt)
    h.fill("sel_e1_dzSignif",dz_signif=sel_vtx.e1.dz/sel_vtx.e1.dzErr,weight=wgt) # double check
    h.fill("sel_e1_chi2",chi2=sel_vtx.e1.trkChi2,weight=wgt)
    h.fill("sel_e1_prob",prob=sel_vtx.e1.trkProb,weight=wgt)
    h.fill("sel_e1_angularRes",angular_res=sel_vtx.e1.angRes,weight=wgt)
    h.fill("sel_e1_trkiso",iso=sel_vtx.e1.trkIso,weight=wgt)
    h.fill("sel_e1_trk_reliso",iso=sel_vtx.e1.trkRelIso,weight=wgt)

    h.fill("sel_e2_dz",dz=sel_vtx.e2.dz,weight=wgt)
    h.fill("sel_e2_dzErr",dz_err=sel_vtx.e2.dzErr,weight=wgt)
    h.fill("sel_e2_dzSignif",dz_signif=sel_vtx.e2.dz/sel_vtx.e2.dzErr,weight=wgt) # double check
    h.fill("sel_e2_chi2",chi2=sel_vtx.e2.trkChi2,weight=wgt)
    h.fill("sel_e2_prob",prob=sel_vtx.e2.trkProb,weight=wgt)
    h.fill("sel_e2_angularRes",angular_res=sel_vtx.e2.angRes,weight=wgt)
    h.fill("sel_e2_trkiso",iso=sel_vtx.e2.trkIso,weight=wgt)
    h.fill("sel_e2_trk_reliso",iso=sel_vtx.e2.trkRelIso,weight=wgt)
    
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
    
    h.fill("sel_vtx_CosThetaColl",cosTheta=sel_vtx.cos_collinear,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV",cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)

    h.fill("sel_vtx_ThetaColl_fromPV_refit",theta=np.degrees(np.arccos(sel_vtx.cos_collinear_fromPV_refit)),weight=wgt)
    h.fill('sel_vtx_ThetaColl_fromPV',theta=np.degrees(np.arccos(sel_vtx.cos_collinear_fromPV)),weight=wgt)


    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_etaErr",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           eta_err=sel_vtx.e1.etaErr,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_etaErr",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           eta_err=sel_vtx.e2.etaErr,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_phiErr",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           phi_err=sel_vtx.e1.phiErr,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_phiErr",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           phi_err=sel_vtx.e2.phiErr,weight=wgt)
    
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_dz",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,dz=sel_vtx.e1.dz,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_dzErr",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dz_err=sel_vtx.e1.dzErr,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_dzSignif",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dz_signif=sel_vtx.e1.dz/sel_vtx.e1.dzErr,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_chi2",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           chi2=sel_vtx.e1.trkChi2,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_prob",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           prob=sel_vtx.e1.trkProb,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_angularRes",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           angular_res=sel_vtx.e1.angRes,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_trkiso",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           iso=sel_vtx.e1.trkIso,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e1_trk_reliso",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           iso=sel_vtx.e1.trkRelIso,weight=wgt)

    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_dz",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dz=sel_vtx.e2.dz,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_dzErr",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dz_err=sel_vtx.e2.dzErr,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_dzSignif",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           dz_signif=sel_vtx.e2.dz/sel_vtx.e2.dzErr,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_chi2",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           chi2=sel_vtx.e2.trkChi2,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_prob",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           prob=sel_vtx.e2.trkProb,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_angularRes",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           angular_res=sel_vtx.e2.angRes,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_trkiso",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           iso=sel_vtx.e2.trkIso,weight=wgt)
    h.fill("sel_vtx_CosThetaColl_fromPV_refit_vs_sel_e2_trk_reliso",cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,\
           iso=sel_vtx.e2.trkRelIso,weight=wgt)

    h.fill('sel_vtx_mass_refit',mass_refit=sel_vtx.refit_m,weight=wgt)
    h.fill('sel_vtx_dR_refit',dR_refit=sel_vtx.refit_dR,weight=wgt)
    