from __future__ import with_statement
import coffea
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema, BaseSchema
from mySchema import MySchema
from coffea import processor
from coffea.nanoevents.methods import vector
from coffea.processor import column_accumulator
import uproot
import awkward as ak
ak.behavior.update(vector.behavior)
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import time
import importlib
import pandas as pd
from XRootD import client
import re
NanoAODSchema.warn_missing_crossrefs = False
import analysisSubroutines as routines
import sys
from collections import defaultdict
from analysisTools import getLumi, deltaPhi

class Skimmer:
    def __init__(self,fileList,max_samples=-1,max_files_per_samp=-1):
        # load in file config
        if type(fileList) == str and ".json" in fileList:
            with open(fileList) as f:
                self.fileList = json.load(f)
        else:
            self.fileList = fileList

        self.sample_names = [] # list of sample names, readable names are generated from data in the fileList json
        self.sample_locs = {} # dictionary mapping sample name to file/directory location
        self.sample_info = {} # dictionary with sample metadata
        self.max_samples = max_samples
        self.max_files_per_samp = max_files_per_samp
        self.nEvents = {}
        self.nEventsProcessed = {}
        self.totalEvents = 0
        self.mode = None
        
        self.loadFiles()
    
    def loadFiles(self):
        loaded = 0
        for sample in self.fileList:
            if self.max_samples > 0 and loaded == self.max_samples:
                break
            mode = sample['type']
            if mode == 'signal':
                name = "sig_{0}_Mchi-{1}_dMchi-{2}_ctau-{3}".format(sample['year'],str(sample['Mchi']).replace(".","p"),str(sample['dMchi']).replace(".","p"),sample['ctau'])
            elif mode == 'bkg':
                name = "bkg_{0}_{1}".format(sample['year'],sample['name'])
            elif mode == 'data':
                name = "data_{0}_{1}".format(sample['year'],sample['name'])
            
            if self.mode is None:
                self.mode = mode
            else:
                if self.mode != mode:
                    print("Error! You're mixing samples of differing types (e.g. signal and bkg, signal and data, etc)")
                    print("Please split up different kinds of samples into different configs")
                    exit()
            
            loc = sample['location']
            if '.root' in loc:
                # if the location is just a single file, load it in
                self.sample_locs[name] = [sample['location']]
            elif 'fileset' in sample.keys():
                self.sample_locs[name] = [f for f in sample['fileset'] if f.split("/")[-1] not in sample['blacklist']]
            else:
                # if the location is a directory, use the xrootd client to get a list of files
                xrdClient = client.FileSystem("root://cmseos.fnal.gov")
                if type(loc) != list:
                    status, flist = xrdClient.dirlist(loc)
                    fullList = ["root://cmsxrootd.fnal.gov/"+loc+"/"+item.name for item in flist if (('.root' in item.name) and (item.name not in sample['blacklist']))]
                else:
                    fullList = []
                    for l in loc:
                        status, flist = xrdClient.dirlist(l)
                        fullList.extend(["root://cmsxrootd.fnal.gov/"+l+"/"+item.name for item in flist if (('.root' in item.name) and (item.name not in sample['blacklist']))])
                if self.max_files_per_samp > 0:
                    fullList = fullList[:self.max_files_per_samp] if len(fullList) > self.max_files_per_samp else fullList
                self.sample_locs[name] = fullList
            
            self.sample_info[name] = sample
            self.sample_names.append(name)
            loaded += 1

    def process(self,treename='ntuples/outT',execr="iterative",workers=4):
        fileset = self.sample_locs
        proc = makeBDTInputs(self.sample_names,self.sample_info,self.sample_locs,mode=self.mode)
        if execr == "iterative":
            executor = processor.IterativeExecutor()
        elif execr == "futures":
            executor = processor.FuturesExecutor(workers=workers)
        else:
            print("Invalid executor type specification!")
            exit
        runner = processor.Runner(executor=executor,schema=MySchema,savemetrics=True)
        accumulator = runner(fileset,
                            treename=treename,
                            processor_instance=proc)
        return accumulator

class makeBDTInputs(processor.ProcessorABC):
    def __init__(self,samples,sampleInfo,fileSet,mode='signal'):
        self.samples = samples
        self.sampleInfo = sampleInfo
        self.sampleLocs = fileSet
        self.mode = mode
    
    def process(self,events):
        samp = events.metadata["dataset"]
        outputs = {}
        info = self.sampleInfo[samp]
        sum_wgt = info["sum_wgt"]
        lumi, unc = getLumi(info['year'])
        xsec = info['xsec']

        # register event weight branch
        events.__setitem__("genWgt",events.genWgt)
        events.__setitem__("sum_wgt",sum_wgt)
        
        events.__setitem__("eventWgt",xsec*lumi*events.genWgt)
        events.__setitem__("eventWgtNorm",xsec*lumi*events.genWgt/sum_wgt)

        # Preselection
        #routines.selectExistingGoodVtx(events)
        routines.selectExistingGoodVtx_noDeltaR(events) # Kyungmin's test
        events.__setitem__("nGoodVtx",ak.count(events.good_vtx.vxy,axis=1))
        events = events[events.nGoodVtx > 0]
        routines.selectBestVertex(events)

        # compute dR and dPhi between selected vertex and jets
        events["sel_vtx","mindRj"] = ak.min(np.sqrt(deltaPhi(events.PFJet.phi,events.sel_vtx.phi)**2 + (events.PFJet.eta-events.sel_vtx.eta)**2),axis=1)
        events["sel_vtx","mindPhiJ"] = ak.min(np.abs(deltaPhi(events.PFJet.phi,events.sel_vtx.phi)),axis=1)

        # computing some signal-only diagnostic quantities
        if info['type'] == "signal":
            e1_match = routines.matchedVertexElectron(events,1)
            e2_match = routines.matchedVertexElectron(events,2)
            genj_phi_pt30 = ak.fill_none(ak.pad_none(events.GenJet.phi[events.GenJet.pt>30],1),999)
            genj_eta_pt30 = ak.fill_none(ak.pad_none(events.GenJet.eta[events.GenJet.pt>30],1),999)
            events["sel_vtx","match"] = ak.where(e1_match*e2_match == -1,2,ak.where(np.abs(e1_match)+np.abs(e2_match) > 0,1,0))
            
            events["GenEle","mindRj"] = ak.min(np.sqrt(deltaPhi(events.PFJet.phi,events.GenEle.phi)**2 + (events.PFJet.eta-events.GenEle.eta)**2),axis=1)
            events["GenEle","mindPhiJ"] = ak.min(np.abs(deltaPhi(events.PFJet.phi,events.GenEle.phi)),axis=1)
            events["GenEle","mindRjGen"] = ak.min(np.sqrt(deltaPhi(genj_phi_pt30,events.GenEle.phi)**2 + (genj_eta_pt30-events.GenEle.eta)**2),axis=1)
            events["GenEle","mindPhiJGen"] = ak.min(ak.where(genj_phi_pt30 != 999,np.abs(deltaPhi(genj_phi_pt30,events.GenEle.phi)),999),axis=1)
            
            events["GenPos","mindRj"] = ak.min(np.sqrt(deltaPhi(events.PFJet.phi,events.GenPos.phi)**2 + (events.PFJet.eta-events.GenPos.eta)**2),axis=1)
            events["GenPos","mindPhiJ"] = ak.min(np.abs(deltaPhi(events.PFJet.phi,events.GenPos.phi)),axis=1)
            events["GenPos","mindRjGen"] = ak.min(np.sqrt(deltaPhi(genj_phi_pt30,events.GenPos.phi)**2 + (genj_eta_pt30-events.GenPos.eta)**2),axis=1)
            events["GenPos","mindPhiJGen"] = ak.min(ak.where(genj_phi_pt30 != 999,np.abs(deltaPhi(genj_phi_pt30,events.GenPos.phi)),999),axis=1)

            events["genEE","mindRj"] = ak.min(np.sqrt(deltaPhi(events.PFJet.phi,events.genEE.phi)**2 + (events.PFJet.eta-events.genEE.eta)**2),axis=1)
            events["genEE","mindPhiJ"] = ak.min(np.abs(deltaPhi(events.PFJet.phi,events.genEE.phi)),axis=1)
            events["genEE","mindRjGen"] = ak.min(np.sqrt(deltaPhi(genj_phi_pt30,events.genEE.phi)**2 + (genj_eta_pt30-events.genEE.eta)**2),axis=1)
            events["genEE","mindPhiJGen"] = ak.min(ak.where(genj_phi_pt30 != 999,np.abs(deltaPhi(genj_phi_pt30,events.genEE.phi)),999),axis=1)

        # filling outputs dictionary
        e1 = events.sel_vtx.e1
        e2 = events.sel_vtx.e2

        outputs['genWgt'] = column_accumulator(events["genWgt"].to_numpy())
        outputs['sum_wgt'] = column_accumulator(events["sum_wgt"].to_numpy())
        
        outputs['wgt'] = column_accumulator(events["eventWgt"].to_numpy())
        outputs['wgt_norm'] = column_accumulator(events["eventWgtNorm"].to_numpy())
        outputs['lead_jet_eta'] = column_accumulator(events.PFJet.eta[:,0].to_numpy())
        outputs['lead_jet_pt'] = column_accumulator(events.PFJet.pt[:,0].to_numpy())
        outputs['jetMETdPhi'] = column_accumulator(np.abs(events.PFJet.METdPhi[:,0]).to_numpy())
        outputs['minJetMETdPhi'] = column_accumulator(ak.min(np.abs(events.PFJet.METdPhi),axis=1).to_numpy())
        outputs['sel_vtx_sign'] = column_accumulator(events.sel_vtx.sign.to_numpy())
        outputs['sel_vtx_chi2'] = column_accumulator(events.sel_vtx.reduced_chi2.to_numpy())
        outputs['sel_vtx_METdPhi'] = column_accumulator(np.abs(events.sel_vtx.METdPhi).to_numpy())
        outputs['sel_vtx_m'] = column_accumulator(events.sel_vtx.m.to_numpy())
        outputs['sel_vtx_dR'] = column_accumulator(events.sel_vtx.dR.to_numpy())
        outputs['sel_vtx_vxy'] = column_accumulator(events.sel_vtx.vxy.to_numpy())
        outputs['sel_vtx_sigmavxy'] = column_accumulator(events.sel_vtx.sigmavxy.to_numpy())
        outputs['sel_vtx_minDxy'] = column_accumulator(np.minimum(np.abs(events.sel_vtx.e1.dxy),np.abs(events.sel_vtx.e2.dxy)).to_numpy())
        outputs['met_leadPt_ratio'] = column_accumulator((events.PFMET.pt/events.PFJet.pt[:,0]).to_numpy())
        outputs["vxy_signif"] = column_accumulator((events.sel_vtx.vxy/events.sel_vtx.sigmavxy).to_numpy())
        outputs["sel_e1_dxy"] = column_accumulator(np.abs(e1.dxy).to_numpy())
        outputs["sel_e1_dxySignif"] = column_accumulator((np.abs(e1.dxy)/e1.dxyErr).to_numpy())
        outputs["sel_e2_dxy"] = column_accumulator(np.abs(e2.dxy).to_numpy())
        outputs["sel_e2_dxySignif"] = column_accumulator((np.abs(e2.dxy)/e2.dxyErr).to_numpy())
        if info['type'] == 'signal':
            outputs['sel_vtx_match'] = column_accumulator(events.sel_vtx.match.to_numpy())
        
            arr_shape = np.shape(events["sum_wgt"].to_numpy())

            ctau = int(samp.split("-")[-1])
            ctau_arr = ctau * np.ones(arr_shape)
            
            outputs["ctau"] = column_accumulator(ctau_arr)

        return outputs

    def postprocess(self, accumulator):
        return accumulator
    