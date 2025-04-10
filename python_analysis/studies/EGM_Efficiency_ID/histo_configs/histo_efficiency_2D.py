import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    # electron
    h.make("gen_ee_vxyCategory_unwgt",'vxyCategories')
    h.make("gen_ee_ptCategory_unwgt",'ptCategories')
    h.make("gen_ee_vxyCategory_vs_gen_ee_ptCategory_unwgt", 'vxyCategories','ptCategories')

    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    if info["type"] == "signal" or info["type"] == "bkg":
        wgt = events.eventWgt/sum_wgt
    else:
        wgt = 1

    h.fill("gen_ee_vxyCategory_unwgt",vxyCat=events.genEE.vxyBin,weight=1)
    h.fill("gen_ee_ptCategory_unwgt",ptCat=events.genEE.ptBin,weight=1)

    h.fill("gen_ee_vxyCategory_vs_gen_ee_ptCategory_unwgt", vxyCat=events.genEE.vxyBin,ptCat=events.genEE.ptBin,weight=1)

