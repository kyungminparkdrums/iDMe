from histobins import *
from hist import Hist
import hist
import numpy as np
import awkward as ak

def make_histograms():
    histograms = {
        # selected vertex displacement
        "sel_vtx_vxy": Hist(samp,cut,vxy,storage=hist.storage.Weight()),
        #"sel_vtx_vxy_zoomed": Hist(samp,cut,vxy_zoomed,storage=hist.storage.Weight()),
        
        # "sel_vtx_sigmavxy": 
        "sel_vtx_minDxy": Hist(samp,cut,ele_dxy,storage=hist.storage.Weight()),
        "sel_vtx_minDz": Hist(samp,cut,ele_dz,storage=hist.storage.Weight()),
        "sel_vtx_deltaDxy": Hist(samp,cut,ele_dxy,storage=hist.storage.Weight()),
        "sel_vtx_deltaDz": Hist(samp,cut,ele_dz,storage=hist.storage.Weight()),
        "sel_vtx_vxySignif": Hist(samp,cut,vtx_vxySignif,storage=hist.storage.Weight()),

        # selected vertex quantities
        "sel_vtx_type": Hist(samp,cut,vtx_type,storage=hist.storage.Weight()),
        "sel_vtx_matchType" : Hist(samp,cut,vtx_matchType,storage=hist.storage.Weight()),
        
        "sel_vtx_chi2": Hist(samp,cut,vtx_chi2,storage=hist.storage.Weight()),
        "sel_vtx_dR": Hist(samp,cut,dR,storage=hist.storage.Weight()),
        "sel_vtx_dEta": Hist(samp,cut,deta,storage=hist.storage.Weight()),
        "sel_vtx_dPhi": Hist(samp,cut,dphi,storage=hist.storage.Weight()),
        "sel_vtx_sign_eta": Hist(samp,cut,sign_eta,storage=hist.storage.Weight()),
        
        "sel_vtx_METdPhi": Hist(samp,cut,dphi,storage=hist.storage.Weight()),
        "sel_vtx_pt": Hist(samp,cut,ele_pt,storage=hist.storage.Weight()),
        "sel_vtx_eta": Hist(samp,cut,ele_eta,storage=hist.storage.Weight()),
        "sel_vtx_phi": Hist(samp,cut,ele_phi,storage=hist.storage.Weight()),
        "sel_vtx_mass": Hist(samp,cut,mass,storage=hist.storage.Weight()),

        # e+/-
        "nLpt_Electron": Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),storage=hist.storage.Weight()),
        "nPF_Electron": Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),storage=hist.storage.Weight()),
        "nElectron": Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),storage=hist.storage.Weight()),

        "nGoodVtx": Hist(samp,cut,hist.axis.Integer(0,10,name="num_good_vtx"),storage=hist.storage.Weight()),
        
        # gen ee dxy, vxy
        

        # 2D selected vertex histos
        "sel_vtx_vxy_vs_mindxy" : Hist(samp,cut,vxy,dxy,storage=hist.storage.Weight()),

        "sel_vtx_vxy_vs_matchType" : Hist(samp,cut,vxy,vtx_matchType,storage=hist.storage.Weight()),

        "sel_vtx_vxy_vs_sel_vtx_chi2" : Hist(samp,cut,vxy,vtx_chi2,storage=hist.storage.Weight()),
        "sel_vtx_vxy_vs_sel_vtx_dR" : Hist(samp,cut,vxy,dR,storage=hist.storage.Weight()),
        "sel_vtx_vxy_vs_sel_vtx_mass" : Hist(samp,cut,vxy,mass,storage=hist.storage.Weight()),
        "sel_vtx_vxy_vs_sel_vtx_METdPhi" : Hist(samp,cut,vxy,dphi,storage=hist.storage.Weight()),
        "sel_vtx_vxy_vs_sel_vtx_vxySignif" : Hist(samp,cut,vxy,vtx_vxySignif,storage=hist.storage.Weight()),
        "sel_vtx_vxy_vs_sel_vtx_pt" : Hist(samp,cut,vxy,ele_pt,storage=hist.storage.Weight()),
        "sel_vtx_vxy_vs_sel_vtx_sign_eta" : Hist(samp,cut,vxy,sign_eta,storage=hist.storage.Weight()),
        
        "sel_vtx_minDxy_vs_matchType" : Hist(samp,cut,ele_dxy,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_minDz_vs_matchType" : Hist(samp,cut,ele_dz,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_deltaDxy_vs_matchType" : Hist(samp,cut,ele_dxy,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_deltaDz_vs_matchType" : Hist(samp,cut,ele_dz,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_vxySignif_vs_matchType" : Hist(samp,cut,vtx_vxySignif,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_chi2_vs_matchType" : Hist(samp,cut,vtx_chi2,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_dR_vs_matchType" : Hist(samp,cut,dR,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_dEta_vs_matchType" : Hist(samp,cut,deta,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_dPhi_vs_matchType" : Hist(samp,cut,dphi,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_METdPhi_vs_matchType" : Hist(samp,cut,dphi,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_pt_vs_matchType" : Hist(samp,cut,ele_pt,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_eta_vs_matchType" : Hist(samp,cut,ele_eta,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_phi_vs_matchType" : Hist(samp,cut,ele_phi,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_mass_vs_matchType" : Hist(samp,cut,mass,vtx_matchType,storage=hist.storage.Weight()),
        "sel_vtx_type_vs_matchType" : Hist(samp,cut,vtx_type,vtx_matchType,storage=hist.storage.Weight()),

        "sel_vtx_sign_eta_vs_matchType" : Hist(samp,cut,sign_eta,vtx_matchType,storage=hist.storage.Weight()),

        "nLpt_Electron_vs_matchType" : Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),vtx_matchType,storage=hist.storage.Weight()),
        "nPF_Electron_vs_matchType" : Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),vtx_matchType,storage=hist.storage.Weight()),
        "nElectron_vs_matchType" : Hist(samp,cut,hist.axis.Integer(0,10,name="num_ele"),vtx_matchType,storage=hist.storage.Weight()),
        "nGoodVtx_vs_matchType" : Hist(samp,cut,hist.axis.Integer(0,10,name="num_good_vtx"),vtx_matchType,storage=hist.storage.Weight()),
        "nGoodVtx_vs_sel_vtx_type" : Hist(samp,cut,hist.axis.Integer(0,10,name="num_good_vtx"),vtx_type,storage=hist.storage.Weight()),
        
        }
    return histograms

subroutines = []

def fillHistos(events,histos,samp,cut,info,sum_wgt=1):
    e1 = events.sel_vtx.e1
    e2 = events.sel_vtx.e2

    min_dxy = np.minimum(np.abs(e1.dxy),np.abs(e2.dxy))
    delta_dxy = np.abs(np.abs(e1.dxy)-np.abs(e2.dxy))

    min_dz = np.minimum(np.abs(e1.dz),np.abs(e2.dz))
    delta_dz = np.abs(np.abs(e1.dz)-np.abs(e2.dz))
    
    max_pfiso = ak.where(e1.PFRelIso<e2.PFRelIso,e2.PFRelIso,e1.PFRelIso)
    
    wgt = events.eventWgt/sum_wgt
    vtx = events.sel_vtx

    histos["sel_vtx_vxy"].fill(samp=samp,cut=cut,vxy=vtx.vxy,weight=wgt)
    #histos["sel_vtx_vxy_zoomed"].fill(samp=samp,cut=cut,vxy_zoomed=vtx.vxy,weight=wgt)
    histos["sel_vtx_minDxy"].fill(samp=samp,cut=cut,dxy=min_dxy,weight=wgt)
    histos["sel_vtx_minDz"].fill(samp=samp,cut=cut,dz=min_dz,weight=wgt)

    histos["sel_vtx_deltaDxy"].fill(samp=samp,cut=cut,dxy=delta_dxy,weight=wgt)
    histos["sel_vtx_deltaDz"].fill(samp=samp,cut=cut,dz=delta_dz,weight=wgt)
    
    histos["sel_vtx_vxySignif"].fill(samp=samp,cut=cut,vxy_signif=vtx.vxy/vtx.sigmavxy,weight=wgt)

    histos["sel_vtx_type"].fill(samp=samp,cut=cut,type=vtx.typ,weight=wgt)

    histos["sel_vtx_chi2"].fill(samp=samp,cut=cut,chi2=vtx.reduced_chi2,weight=wgt)
    histos["sel_vtx_dR"].fill(samp=samp,cut=cut,dr=vtx.dR,weight=wgt)

    histos["sel_vtx_dEta"].fill(samp=samp,cut=cut,deta=np.abs(e1.eta-e2.eta),weight=wgt)
    histos["sel_vtx_dPhi"].fill(samp=samp,cut=cut,dphi=np.abs(e1.phi-e2.phi),weight=wgt)

    histos["sel_vtx_sign_eta"].fill(samp=samp,cut=cut,sign_eta=(e1.eta*e2.eta)/np.abs(e1.eta*e2.eta),weight=wgt)
    
    histos["sel_vtx_METdPhi"].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),weight=wgt)

    histos["sel_vtx_pt"].fill(samp=samp,cut=cut,pt=vtx.pt,weight=wgt)
    histos["sel_vtx_eta"].fill(samp=samp,cut=cut,eta=vtx.eta,weight=wgt)
    histos["sel_vtx_phi"].fill(samp=samp,cut=cut,phi=vtx.phi,weight=wgt)
    histos["sel_vtx_mass"].fill(samp=samp,cut=cut,mass=vtx.m,weight=wgt)

    histos['nLpt_Electron'].fill(samp=samp,cut=cut,num_ele=ak.count(events.LptElectron.pt,axis=1),weight=wgt)
    histos['nPF_Electron'].fill(samp=samp,cut=cut,num_ele=ak.count(events.Electron.pt,axis=1),weight=wgt)
    histos['nElectron'].fill(samp=samp,cut=cut,num_ele=ak.count(events.LptElectron.pt,axis=1)+ak.count(events.Electron.pt,axis=1),weight=wgt)

    histos['nGoodVtx'].fill(samp=samp,cut=cut,num_good_vtx=ak.count(events.good_vtx.vxy,axis=1),weight=wgt)
    
    histos['sel_vtx_vxy_vs_mindxy'].fill(samp=samp,cut=cut,vxy=vtx.vxy,dxy=min_dxy,weight=wgt)

    histos['sel_vtx_vxy_vs_sel_vtx_chi2'].fill(samp=samp,cut=cut,vxy=vtx.vxy,chi2=vtx.reduced_chi2,weight=wgt)
    histos['sel_vtx_vxy_vs_sel_vtx_dR'].fill(samp=samp,cut=cut,vxy=vtx.vxy,dr=vtx.dR,weight=wgt)
    histos['sel_vtx_vxy_vs_sel_vtx_mass'].fill(samp=samp,cut=cut,vxy=vtx.vxy,mass=vtx.m,weight=wgt)
    histos['sel_vtx_vxy_vs_sel_vtx_METdPhi'].fill(samp=samp,cut=cut,vxy=vtx.vxy,dphi=np.abs(vtx.METdPhi),weight=wgt)
    histos['sel_vtx_vxy_vs_sel_vtx_vxySignif'].fill(samp=samp,cut=cut,vxy=vtx.vxy,vxy_signif=vtx.vxy/vtx.sigmavxy,weight=wgt)
    histos['sel_vtx_vxy_vs_sel_vtx_pt'].fill(samp=samp,cut=cut,vxy=vtx.vxy,pt=vtx.pt,weight=wgt)
    histos['sel_vtx_vxy_vs_sel_vtx_sign_eta'].fill(samp=samp,cut=cut,vxy=vtx.vxy,sign_eta=(e1.eta*e2.eta)/np.abs(e1.eta*e2.eta),weight=wgt)
    
    if info["type"] == "signal":
        histos["sel_vtx_matchType"].fill(samp=samp,cut=cut,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_vxy_vs_matchType'].fill(samp=samp,cut=cut,vxy=vtx.vxy,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_minDxy_vs_matchType'].fill(samp=samp,cut=cut,dxy=min_dxy,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_minDz_vs_matchType'].fill(samp=samp,cut=cut,dz=min_dz,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_deltaDxy_vs_matchType'].fill(samp=samp,cut=cut,dxy=delta_dxy,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_deltaDz_vs_matchType'].fill(samp=samp,cut=cut,dz=delta_dz,mtype=vtx.match,weight=wgt)

        histos['sel_vtx_vxySignif_vs_matchType'].fill(samp=samp,cut=cut,vxy_signif=vtx.vxy/vtx.sigmavxy,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_chi2_vs_matchType'].fill(samp=samp,cut=cut,chi2=vtx.reduced_chi2,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_dR_vs_matchType'].fill(samp=samp,cut=cut,dr=vtx.dR,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_dEta_vs_matchType'].fill(samp=samp,cut=cut,deta=np.abs(e1.eta-e2.eta),mtype=vtx.match,weight=wgt)
        histos['sel_vtx_dPhi_vs_matchType'].fill(samp=samp,cut=cut,dphi=np.abs(e1.phi-e2.phi),mtype=vtx.match,weight=wgt)

        histos["sel_vtx_sign_eta_vs_matchType"].fill(samp=samp,cut=cut,sign_eta=(e1.eta*e2.eta)/np.abs(e1.eta*e2.eta),mtype=vtx.match,weight=wgt)
        
        histos['sel_vtx_METdPhi_vs_matchType'].fill(samp=samp,cut=cut,dphi=np.abs(vtx.METdPhi),mtype=vtx.match,weight=wgt)
        histos['sel_vtx_pt_vs_matchType'].fill(samp=samp,cut=cut,pt=vtx.pt,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_eta_vs_matchType'].fill(samp=samp,cut=cut,eta=vtx.eta,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_phi_vs_matchType'].fill(samp=samp,cut=cut,phi=vtx.phi,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_mass_vs_matchType'].fill(samp=samp,cut=cut,mass=vtx.m,mtype=vtx.match,weight=wgt)
        histos['sel_vtx_type_vs_matchType'].fill(samp=samp,cut=cut,type=vtx.typ,mtype=vtx.match,weight=wgt)

        histos['nLpt_Electron_vs_matchType'].fill(samp=samp,cut=cut,num_ele=ak.count(events.LptElectron.pt,axis=1),mtype=vtx.match,weight=wgt)
        histos['nPF_Electron_vs_matchType'].fill(samp=samp,cut=cut,num_ele=ak.count(events.Electron.pt,axis=1),mtype=vtx.match,weight=wgt)
        histos['nElectron_vs_matchType'].fill(samp=samp,cut=cut,num_ele=ak.count(events.LptElectron.pt,axis=1)+ak.count(events.Electron.pt,axis=1)\
                                              ,mtype=vtx.match,weight=wgt)
        histos['nGoodVtx_vs_matchType'].fill(samp=samp,cut=cut,num_good_vtx=ak.count(events.good_vtx.vxy,axis=1),mtype=vtx.match,weight=wgt)
        histos['nGoodVtx_vs_sel_vtx_type'].fill(samp=samp,cut=cut,num_good_vtx=ak.count(events.good_vtx.vxy,axis=1),type=vtx.typ,weight=wgt)
