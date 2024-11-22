import numpy as np
import awkward as ak
from myHisto import myHisto

# General Purpose
#samp = StrCategory([],name="samp",label="Sample Name",growth=True)
#cut = StrCategory([],name="cut",label="Cut Applied",growth=True)

# functions to make histograms
"""def parse_axis(a):
    name = a[0]
    if type(a[1]) == list:
        assert len(a) == 2
        if type(a[1][0]) == str:
            axis = StrCategory(a[1],name=name,label=name)
        else:
            axis = IntCategory(a[1],name=name,label=name)
    else:
        assert len(a) == 4
        axis = Regular(a[1],a[2],a[3],name=name,label=name)
    return axis

def histo(name,*args):
    axes = [samp,cut]
    for ax in args:
        axes.append(parse_axis(ax))
    return Hist(*axes,storage=hist.storage.Weight())"""
    
def make_histograms(info):
    h = myHisto()

    h.make("gen_met",'met')
    #histograms = {
    #    # quantities associated w/ gen objects
    #    "gen_met" : histo("gen_met",('met',100,50,300)),
    #    "gen_met_noWgt" : histo("gen_met",('met',100,50,300)),
    #    "gen_leadjet_pt_noWgt" : histo("gen_leadjet_pt_noWgt", ('pt',100,50,300)),
    #    "gen_leadjet_pt" : histo("gen_leadjet_pt", ('pt',100,50,300)),
    #    "gen_met_vs_leadjet_pt_noWgt" : histo("gen_met_vs_leadjet_pt_noWgt",('met',100,50,300),('pt',100,50,300)),
    #    "gen_met_vs_leadjet_pt" : histo("gen_met_vs_leadjet_pt",('met',100,50,300),('pt',100,50,300))
    #}
    #return histograms
    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    wgt = events.eventWgt/sum_wgt
    h.samp = samp
    h.cut = cut
    
    h.fill("gen_met",met=events.GenMET.pt,weight=1)
    
    #if info["type"] == "signal":
        #histos["gen_met_noWgt"].fill(samp=samp,cut=cut,met=events.GenMET.pt,weight=1)
    #    histos["gen_met"].fill(samp=samp,cut=cut,met=events.GenMET.pt,weight=wgt)
    #    histos["gen_leadjet_pt_noWgt"].fill(samp=samp,cut=cut,pt=events.GenJet.pt[:,0],weight=1)
    #    histos["gen_leadjet_pt"].fill(samp=samp,cut=cut,pt=events.GenJet.pt[:,0],weight=wgt)
    #    histos["gen_met_vs_leadjet_pt_noWgt"].fill(samp=samp,cut=cut,met=events.GenMET.pt,pt=events.GenJet.pt[:,0],weight=1)
    #    histos["gen_met_vs_leadjet_pt"].fill(samp=samp,cut=cut,met=events.GenMET.pt,pt=events.GenJet.pt[:,0],weight=wgt)