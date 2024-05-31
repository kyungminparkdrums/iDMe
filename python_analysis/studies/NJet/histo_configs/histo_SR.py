from histobins import *
from hist import Hist
import hist
import numpy as np
import awkward as ak

def make_histograms():
    histograms = {
        #"bdt_score" : Hist(samp,cut,bdtScore,storage=hist.storage.Weight()),
        
        # selected vertex displacement
        "sel_vtx_vxy": Hist(samp,cut,vxy,storage=hist.storage.Weight()),
        "sel_vtx_vxy_projected": Hist(samp,cut,vxy_projected,storage=hist.storage.Weight()),
        "sel_vtx_minDxy": Hist(samp,cut,ele_dxy,storage=hist.storage.Weight()),
        "sel_vtx_vxySignif": Hist(samp,cut,vtx_vxySignif,storage=hist.storage.Weight()),

        "cos_collinear": Hist(samp,cut,cos_collinear,storage=hist.storage.Weight()),   

        # selected vertex quantities
        "sel_vtx_type": Hist(samp,cut,vtx_type,storage=hist.storage.Weight()),
        "sel_vtx_chi2": Hist(samp,cut,vtx_chi2,storage=hist.storage.Weight()),
        "sel_vtx_dR": Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "sel_vtx_sign_eta": Hist(samp,cut,sign_eta,storage=hist.storage.Weight()),
        "sel_vtx_prod_eta": Hist(samp,cut,prod_eta,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi": Hist(samp,cut,dphi,storage=hist.storage.Weight()),
        "sel_vtx_pt": Hist(samp,cut,ele_pt,storage=hist.storage.Weight()),
        "sel_vtx_mass": Hist(samp,cut,mass,storage=hist.storage.Weight()),

        "sel_vtx_pt_e1_over_pt_e2": Hist(samp,cut,ratio,storage=hist.storage.Weight()),
        "delta_dxy_over_maxdxy": Hist(samp,cut,ratio,storage=hist.storage.Weight()),

        "sel_vtx_pt_over_m": Hist(samp,cut,sel_vtx_pt_over_m,storage=hist.storage.Weight()),
        "sel_vtx_log_dEta_over_dPhi": Hist(samp,cut,ratio_big,storage=hist.storage.Weight()),

        # e+/-
        "nLpt_Electron": Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),storage=hist.storage.Weight()),
        "nPF_Electron": Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),storage=hist.storage.Weight()),
        "nElectron": Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),storage=hist.storage.Weight()),

        "nGoodVtx": Hist(samp,cut,hist.axis.Integer(0,10,name="num_good_vtx"),storage=hist.storage.Weight()),

        # ISR Jets
        "met_over_lead_jet_pt": Hist(samp,cut,met_over_pt,storage=hist.storage.Weight()),
        "lead_jet_met_dPhi" : Hist(samp,cut,dphi,storage=hist.storage.Weight()),
        "min_jet_met_dPhi" : Hist(samp,cut,dphi,storage=hist.storage.Weight()),
        "MET_pt" : Hist(samp,cut,met_pt,storage=hist.storage.Weight()),
        "nJets" : Hist(samp,cut,njets,storage=hist.storage.Weight()),
        }
    return histograms

subroutines = []
    

def fillHistos(events,histos,samp,cut,info,sum_wgt=1):
    e1 = events.sel_vtx.e1
    e2 = events.sel_vtx.e2

    min_dxy = np.minimum(np.abs(e1.dxy),np.abs(e2.dxy))
    max_dxy = np.maximum(np.abs(e1.dxy),np.abs(e2.dxy))
    delta_dxy = np.abs(np.abs(e1.dxy)-np.abs(e2.dxy))

    min_dz = np.minimum(np.abs(e1.dz),np.abs(e2.dz))
    delta_dz = np.abs(np.abs(e1.dz)-np.abs(e2.dz))
    
    max_pfiso = ak.where(e1.PFRelIso<e2.PFRelIso,e2.PFRelIso,e1.PFRelIso)

    if info["type"] == "signal" or info["type"] == "bkg":
        wgt = events.eventWgt/sum_wgt
    else:
        wgt = 1

    vtx = events.sel_vtx

    #histos["bdt_score"].fill(samp=samp,cut=cut,score=events.BDTScore,weight=wgt)
    
    histos["sel_vtx_vxy"].fill(samp=samp,cut=cut,vxy=vtx.vxy,weight=wgt)
    histos["sel_vtx_vxy_projected"].fill(samp=samp,cut=cut,vxy_projected=events.sel_vtx.projectedLxy,weight=wgt)
    histos["sel_vtx_minDxy"].fill(samp=samp,cut=cut,dxy=min_dxy,weight=wgt)
    histos["sel_vtx_vxySignif"].fill(samp=samp,cut=cut,vxy_signif=vtx.vxy/vtx.sigmavxy,weight=wgt)

    histos['cos_collinear'].fill(samp=samp,cut=cut,cos_collinear=events.sel_vtx.cos_collinear,weight=wgt)

    histos["sel_vtx_type"].fill(samp=samp,cut=cut,type=vtx.typ,weight=wgt)
    histos["sel_vtx_chi2"].fill(samp=samp,cut=cut,chi2=vtx.reduced_chi2,weight=wgt)
    histos["sel_vtx_dR"].fill(samp=samp,cut=cut,dr=vtx.dR,weight=wgt)
    histos["sel_vtx_sign_eta"].fill(samp=samp,cut=cut,sign_eta=(e1.eta*e2.eta)/np.abs(e1.eta*e2.eta),weight=wgt)
    histos["sel_vtx_prod_eta"].fill(samp=samp,cut=cut,prod_eta=e1.eta*e2.eta,weight=wgt)
    histos["sel_vtx_METdPhi"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),weight=wgt)
    histos["sel_vtx_pt"].fill(samp=samp,cut=cut,pt=vtx.pt,weight=wgt)
    histos["sel_vtx_mass"].fill(samp=samp,cut=cut,mass=vtx.m,weight=wgt)

    histos["sel_vtx_pt_e1_over_pt_e2"].fill(samp=samp,cut=cut,ratio=np.minimum(e1.pt, e2.pt)/np.maximum(e1.pt, e2.pt),weight=wgt)
    histos["delta_dxy_over_maxdxy"].fill(samp=samp,cut=cut,ratio=delta_dxy/max_dxy,weight=wgt)

    histos["sel_vtx_pt_over_m"].fill(samp=samp,cut=cut,sel_vtx_pt_over_m=vtx.pt/vtx.m,weight=wgt)
    histos["sel_vtx_log_dEta_over_dPhi"].fill(samp=samp,cut=cut,ratio_big=np.log10(np.abs(e1.eta-e2.eta)/np.abs(e1.phi-e2.phi)),weight=wgt)

    histos['nLpt_Electron'].fill(samp=samp,cut=cut,num_ele=ak.count(events.LptElectron.pt,axis=1),weight=wgt)
    histos['nPF_Electron'].fill(samp=samp,cut=cut,num_ele=ak.count(events.Electron.pt,axis=1),weight=wgt)
    histos['nElectron'].fill(samp=samp,cut=cut,num_ele=ak.count(events.LptElectron.pt,axis=1)+ak.count(events.Electron.pt,axis=1),weight=wgt)

    histos['nGoodVtx'].fill(samp=samp,cut=cut,num_good_vtx=ak.count(events.good_vtx.vxy,axis=1),weight=wgt)

    histos["met_over_lead_jet_pt"].fill(samp=samp,cut=cut,met_over_pt=events.PFMET.pt/events.PFJet.pt[:,0],weight=wgt)
    histos["lead_jet_met_dPhi"].fill(samp=samp,cut=cut,dphi=np.abs(events.PFJet.METdPhi[:,0]),weight=wgt)
    histos["min_jet_met_dPhi"].fill(samp=samp,cut=cut,dphi=ak.min(np.abs(events.PFJet.METdPhi),axis=1),weight=wgt)
    histos["MET_pt"].fill(samp=samp,cut=cut,met_pt=events.PFMET.pt,weight=wgt)
    histos["nJets"].fill(samp=samp,cut=cut,njets=ak.count(events.PFJet.pt,axis=1),weight=wgt)