import numpy as np
import awkward as ak
from myHisto import myHisto

def make_histograms(info):
    h = myHisto()
    
    # matched reco electron, reco variables
    h.make("PFMET",'met')
    h.make("nJets",'nJets')
    h.make("nBjets_loose",'nJets')
    h.make("nBjets_med",'nJets')
    h.make("nBjets_tight",'nJets')
    
    # gen info
    if info['type'] == 'signal':
        h.make("gen_dR",'dR')
        h.make("gen_MET",'met')
        h.make("gen_vxy1",'vxy1')
        h.make("gen_vxy10",'vxy10')
        h.make("gen_vxy100",'vxy100')
        h.make("gen_leadpT",'ele_pt_fine')
        h.make("gen_eeMETdPhi",'abs_dphi')
        h.make("gen_jetMETdPhi",'abs_dphi')
        h.make("gen_jetMETratio",'jetMETratio')
        h.make('gen_lead_jet_pt','jet_pt')
        h.make('gen_lead_jet_eta','eta')
        
        h.make("match_ele_dxy",'dxy')

    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    if info["type"] == "signal" or info["type"] == "bkg":
        wgt = events.eventWgt/sum_wgt
    else:
        wgt = 1

    ### FILLING HISTOGRAMS ###
    h.fill("PFMET",met=events.PFMET.pt,weight=wgt)
    h.fill("nJets",nJets=ak.count(events.PFJet.pt,axis=1),weight=wgt)
    h.fill("nBjets_loose",nJets=ak.count_nonzero(events.PFJet.passLooseID,axis=1),weight=wgt)
    h.fill("nBjets_med",nJets=ak.count_nonzero(events.PFJet.passMedID,axis=1),weight=wgt)
    h.fill("nBjets_tight",nJets=ak.count_nonzero(events.PFJet.passTightID,axis=1),weight=wgt)
    #
    if info['type'] == 'signal':
        h.fill("gen_dR",dR=events.genEE.dr,weight=wgt)
        h.fill("gen_MET",met=events.GenMET.pt,weight=wgt)
        h.fill("gen_vxy1",vxy=events.genEE.vxy,weight=wgt)
        h.fill("gen_vxy10",vxy=events.genEE.vxy,weight=wgt)
        h.fill("gen_vxy100",vxy=events.genEE.vxy,weight=wgt)
        h.fill("gen_leadpT",pt=np.maximum(events.GenEle.pt,events.GenPos.pt),weight=wgt)
        h.fill("gen_eeMETdPhi",abs_dphi=np.abs(events.genEE.METdPhi),weight=wgt)
        h.fill("gen_jetMETdPhi",abs_dphi=np.abs(events.GenJet.METdPhi[:,0]),weight=wgt)
        h.fill("gen_jetMETratio",jetMETratio=events.GenJet.pt[:,0]/events.GenMET.pt,weight=wgt)
        h.fill('gen_lead_jet_pt',jet_pt=events.GenJet.pt[:,0],weight=wgt)
        h.fill('gen_lead_jet_eta',eta=events.GenJet.eta[:,0],weight=wgt)
        
        h.fill("match_ele_dxy",dxy=np.abs(ak.flatten(events.Electron[events.Electron.genMatched].dxy)),weight=1)
        h.fill("match_ele_dxy",dxy=np.abs(ak.flatten(events.LptElectron[events.LptElectron.genMatched].dxy)),weight=1)
               
        