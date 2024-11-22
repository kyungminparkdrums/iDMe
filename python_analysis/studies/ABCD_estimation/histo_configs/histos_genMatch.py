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
    
    ## ##
    
    h.make("sel_vtx_vx_vs_vy","vx","vy")
    h.make("PVx_vs_PVy","pvx","pvy")
    
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

    # Test "fake" gen-matching and vtx sign
    h.make('sel_vtx_sign_vs_sel_vtx_type','vtx_sign','vtx_type')
    h.make('sel_vtx_sign_vs_sel_e1_pt_gen_res','vtx_sign','ele_ptRes')
    h.make('sel_vtx_sign_vs_sel_e2_pt_gen_res','vtx_sign','ele_ptRes')
    h.make('sel_vtx_sign_vs_sel_e1_gen_dR','vtx_sign','dR_zoom')
    h.make('sel_vtx_sign_vs_sel_e2_gen_dR','vtx_sign','dR_zoom')

    h.make('sel_vtx_sign_vs_sel_vtx_thetaColl_fromPV_refit_rad','vtx_sign','theta_rad')

    h.make('sel_vtx_sign_vs_sel_vtx_mindxy_refit','vtx_sign','dxy')
    h.make('sel_vtx_sign_vs_sel_vtx_vxy10','vtx_sign','vxy10')
    h.make('sel_vtx_sign_vs_sel_vtx_vxy1','vtx_sign','vxy1')

    h.make('sel_vtx_dR_vs_sel_vtx_dR_gen','dR','dR_gen')
    h.make('sel_vtx_sign_vs_sel_vtx_pt','vtx_sign','ele_pt')
    h.make('e1_sign_reco_gen','ele_sign','gen_ele_sign')
    h.make('e2_sign_reco_gen','ele_sign','gen_ele_sign')

    h.make('gen_e_pT','ele_pt')
    h.make('gen_p_pT','ele_pt')

    h.make('sel_e1_pt_gen_res','ele_ptRes')
    h.make('sel_e2_pt_gen_res','ele_ptRes')

    h.make('sel_e1_gen_dR','dR_zoom')
    h.make('sel_e2_gen_dR','dR_zoom')

    h.make('e1_sign_prod_vs_e1_pt','ele_sign','ele_pt')
    h.make('e1_sign_prod_vs_e1_eta','ele_sign','eta')
    h.make('e1_sign_prod_vs_e1_phi','ele_sign','phi')
    h.make('e1_sign_prod_vs_e1_refit_dxy','ele_sign','dxy')

    h.make('e2_sign_prod_vs_e2_pt','ele_sign','ele_pt')
    h.make('e2_sign_prod_vs_e2_eta','ele_sign','eta')
    h.make('e2_sign_prod_vs_e2_phi','ele_sign','phi')
    h.make('e2_sign_prod_vs_e2_refit_dxy','ele_sign','dxy')

    
    return h

subroutines = []

def dxy_custom(vx, vy, beamspot_x, beamspot_y, px, py, pt):
    dxy = ((-(vx - beamspot_x) * py) + (vy - beamspot_y)*px)/pt

    return dxy

def get_dR(gen_eta, gen_phi, reco_eta, reco_phi):
    dR_sq = (gen_eta - reco_eta)**2 + (gen_phi - reco_phi)**2
    return np.sqrt(dR_sq)
    

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
        ####
        dR_e1_gen_e = get_dR(events.GenEle.eta, events.GenEle.phi, sel_vtx.e1.eta, sel_vtx.e1.phi)
        dR_e1_gen_p = get_dR(events.GenPos.eta, events.GenPos.phi, sel_vtx.e1.eta, sel_vtx.e1.phi)
        dR_e2_gen_e = get_dR(events.GenEle.eta, events.GenEle.phi, sel_vtx.e2.eta, sel_vtx.e2.phi)
        dR_e2_gen_p = get_dR(events.GenPos.eta, events.GenPos.phi, sel_vtx.e2.eta, sel_vtx.e2.phi)

        mask_e1_e = dR_e1_gen_e < dR_e1_gen_p
        dR_e1_gen = ak.fill_none(ak.mask(dR_e1_gen_e, mask_e1_e),0)+ak.fill_none(ak.mask(dR_e1_gen_p, ~mask_e1_e),0)

        e1_genEle_pt_res = np.abs(sel_vtx.e1.pt - events.GenEle.pt)/events.GenEle.pt
        e1_genPos_pt_res = np.abs(sel_vtx.e1.pt - events.GenPos.pt)/events.GenPos.pt
        
        sel_e1_pt_gen_res = ak.fill_none(ak.mask(e1_genEle_pt_res, mask_e1_e),0)+ak.fill_none(ak.mask(e1_genPos_pt_res, ~mask_e1_e),0)
        gen_e1_sign =  ak.fill_none(ak.mask((-1)*np.ones(len(events)), mask_e1_e),0)+ak.fill_none(ak.mask(np.ones(len(events)), ~mask_e1_e),0)
        
        mask_e2_e = dR_e2_gen_e < dR_e2_gen_p
        dR_e2_gen = ak.fill_none(ak.mask(dR_e2_gen_e, mask_e2_e),0)+ak.fill_none(ak.mask(dR_e2_gen_p, ~mask_e2_e),0)

        e2_genEle_pt_res = np.abs(sel_vtx.e2.pt - events.GenEle.pt)/events.GenEle.pt
        e2_genPos_pt_res = np.abs(sel_vtx.e2.pt - events.GenPos.pt)/events.GenPos.pt
        
        sel_e2_pt_gen_res = ak.fill_none(ak.mask(e2_genEle_pt_res, mask_e2_e),0)+ak.fill_none(ak.mask(e2_genPos_pt_res, ~mask_e2_e),0)
        gen_e2_sign =  ak.fill_none(ak.mask((-1)*np.ones(len(events)), mask_e2_e),0)+ak.fill_none(ak.mask(np.ones(len(events)), ~mask_e2_e),0)

        e1_sign_prod = sel_vtx.e1.charge * gen_e1_sign
        e2_sign_prod = sel_vtx.e2.charge * gen_e2_sign

        h.fill('e1_sign_prod_vs_e1_pt',sign=e1_sign_prod,pt=sel_vtx.e1.pt,weight=wgt)
        h.fill('e1_sign_prod_vs_e1_eta',sign=e1_sign_prod,eta=sel_vtx.e1.eta,weight=wgt)
        h.fill('e1_sign_prod_vs_e1_phi',sign=e1_sign_prod,phi=sel_vtx.e1.phi,weight=wgt)
        h.fill('e1_sign_prod_vs_e1_refit_dxy',sign=e1_sign_prod,dxy=sel_vtx.e1.refit_dxy,weight=wgt)

        h.fill('e2_sign_prod_vs_e2_pt',sign=e2_sign_prod,pt=sel_vtx.e2.pt,weight=wgt)
        h.fill('e2_sign_prod_vs_e2_eta',sign=e2_sign_prod,eta=sel_vtx.e2.eta,weight=wgt)
        h.fill('e2_sign_prod_vs_e2_phi',sign=e2_sign_prod,phi=sel_vtx.e2.phi,weight=wgt)
        h.fill('e2_sign_prod_vs_e2_refit_dxy',sign=e2_sign_prod,dxy=sel_vtx.e2.refit_dxy,weight=wgt)
        
        h.fill('sel_vtx_sign_vs_sel_vtx_type',sign=sel_vtx.sign,vtype=sel_vtx.typ,weight=wgt)
        h.fill('sel_vtx_sign_vs_sel_e1_pt_gen_res',sign=sel_vtx.sign,ptres=sel_e1_pt_gen_res,weight=wgt)
        h.fill('sel_vtx_sign_vs_sel_e2_pt_gen_res',sign=sel_vtx.sign,ptres=sel_e2_pt_gen_res,weight=wgt)
        h.fill('sel_vtx_sign_vs_sel_e1_gen_dR',sign=sel_vtx.sign,dR=dR_e1_gen,weight=wgt)
        h.fill('sel_vtx_sign_vs_sel_e2_gen_dR',sign=sel_vtx.sign,dR=dR_e2_gen,weight=wgt)
        h.fill('sel_vtx_sign_vs_sel_vtx_thetaColl_fromPV_refit_rad',sign=sel_vtx.sign,theta_rad=np.arccos(sel_vtx.cos_collinear_fromPV),weight=wgt)

        h.fill('sel_vtx_sign_vs_sel_vtx_mindxy_refit',sign=sel_vtx.sign,dxy=np.minimum(sel_vtx.e1.refit_dxy, sel_vtx.e2.refit_dxy),weight=wgt)
        h.fill('sel_vtx_sign_vs_sel_vtx_vxy1',sign=sel_vtx.sign,vxy=sel_vtx.vxy,weight=wgt)
        h.fill('sel_vtx_sign_vs_sel_vtx_vxy10',sign=sel_vtx.sign,vxy=sel_vtx.vxy,weight=wgt)

        h.fill('sel_vtx_dR_vs_sel_vtx_dR_gen',dR=sel_vtx.dR,dR_gen=get_dR(events.GenEle.eta, events.GenEle.phi, events.GenPos.eta, events.GenPos.phi),weight=wgt)
        h.fill('sel_vtx_sign_vs_sel_vtx_pt',sign=sel_vtx.sign,pt=sel_vtx.pt,weight=wgt)
        h.fill('e1_sign_reco_gen',sign=sel_vtx.e1.charge,gen_sign=gen_e1_sign,weight=wgt)
        h.fill('e2_sign_reco_gen',sign=sel_vtx.e2.charge,gen_sign=gen_e2_sign,weight=wgt)

        h.fill('gen_e_pT',pt=events.GenEle.pt,weight=wgt)
        h.fill('gen_p_pT',pt=events.GenPos.pt,weight=wgt)

        h.fill('sel_e1_pt_gen_res',ptres=sel_e1_pt_gen_res,weight=wgt)
        h.fill('sel_e2_pt_gen_res',ptres=sel_e2_pt_gen_res,weight=wgt)

        h.fill('sel_e1_gen_dR',dR=dR_e1_gen,weight=wgt)
        h.fill('sel_e2_gen_dR',dR=dR_e2_gen,weight=wgt)

        ####
        
        h.fill("sel_vtx_GenCosThetaColl",cosTheta=sel_vtx.gen_cos_collinear_fromPV,weight=wgt)
        h.fill("sel_vtx_GenThetaColl",gen_theta=np.degrees(np.arccos(sel_vtx.gen_cos_collinear_fromPV)),weight=wgt)
        h.fill("sel_vtx_GenCosThetaColl_vs_sel_vtx_CosThetaColl_fromPV",cosTheta=sel_vtx.gen_cos_collinear_fromPV,\
               cosTheta_fromPV=sel_vtx.cos_collinear_fromPV,weight=wgt)
        h.fill("sel_vtx_GenCosThetaColl_vs_sel_vtx_CosThetaColl_fromPV_refit",cosTheta=sel_vtx.gen_cos_collinear_fromPV,\
               cosTheta_fromPV_refit=sel_vtx.cos_collinear_fromPV_refit,weight=wgt)
        h.fill("sel_vtx_GenThetaColl_vs_sel_vtx_ThetaColl_fromPV",gen_theta=np.degrees(np.arccos(sel_vtx.gen_cos_collinear_fromPV)),\
               theta=np.degrees(np.arccos(sel_vtx.cos_collinear_fromPV)),weight=wgt)
        h.fill("sel_vtx_GenThetaColl_vs_sel_vtx_ThetaColl_fromPV_refit",gen_theta=np.degrees(np.arccos(sel_vtx.gen_cos_collinear_fromPV)),\
               theta=np.degrees(np.arccos(sel_vtx.cos_collinear_fromPV_refit)),weight=wgt)
        h.fill("sel_vtx_GenThetaColl_rad_vs_sel_vtx_ThetaColl_fromPV_refit_rad",gen_theta_rad=np.arccos(sel_vtx.gen_cos_collinear_fromPV),\
               theta_rad=np.arccos(sel_vtx.cos_collinear_fromPV_refit),weight=wgt)