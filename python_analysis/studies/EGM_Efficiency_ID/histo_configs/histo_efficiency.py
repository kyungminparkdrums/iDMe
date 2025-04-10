import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    # electron
    h.make("gen_ele_vxy20",'vxy20')
    h.make("gen_ele_vxy100",'vxy100')
    h.make("gen_ele_pt",'ele_pt')

    h.make("gen_ele_vxy20_vs_reco_type",'vxy20', 'reco_type')
    h.make("gen_ele_vxy100_vs_reco_type",'vxy100', 'reco_type')
    h.make("gen_ele_pt_vs_reco_type",'ele_pt', 'reco_type')

    # pos
    h.make("gen_pos_vxy20",'vxy20')
    h.make("gen_pos_vxy100",'vxy100')
    h.make("gen_pos_pt",'ele_pt')

    h.make("gen_pos_vxy20_vs_reco_type",'vxy20', 'reco_type')
    h.make("gen_pos_vxy100_vs_reco_type",'vxy100', 'reco_type')
    h.make("gen_pos_pt_vs_reco_type",'ele_pt', 'reco_type')

    # ele
    #h.make("gen_vxyCategory_unwgt",'vxyCategories')
    #h.make("gen_ptCategory_unwgt",'ptCategories')

    # ee
    h.make("gen_ee_vxy20_vs_reco_type",'vxy20', 'vtx_reco_type')
    h.make("gen_ee_vxy100_vs_reco_type",'vxy100', 'vtx_reco_type')
    h.make("gen_ee_pt_vs_reco_type",'vtx_pt', 'vtx_reco_type')
    
    return h

subroutines = []

def fillHistos(events,h,samp,cut,info,sum_wgt=1):
    h.samp = samp
    h.cut = cut
    if info["type"] == "signal" or info["type"] == "bkg":
        wgt = events.eventWgt/sum_wgt
    else:
        wgt = 1

    #sel_vtx = events.sel_vtx

    h.fill("gen_ele_vxy20",vxy=events.GenEle.vxy,weight=1)
    h.fill("gen_ele_vxy100",vxy=events.GenEle.vxy,weight=1)
    h.fill("gen_ele_pt",pt=events.GenEle.pt,weight=1)

    h.fill("gen_ele_vxy20_vs_reco_type", vxy=events.GenEle.vxy, reco_type=events.GenEle.matchType, weight=1)
    h.fill("gen_ele_vxy100_vs_reco_type", vxy=events.GenEle.vxy, reco_type=events.GenEle.matchType, weight=1)
    h.fill("gen_ele_pt_vs_reco_type", pt=events.GenEle.pt, reco_type=events.GenEle.matchType, weight=1)

    h.fill("gen_pos_vxy20",vxy=events.GenPos.vxy,weight=1)
    h.fill("gen_pos_vxy100",vxy=events.GenPos.vxy,weight=1)
    h.fill("gen_pos_pt",pt=events.GenPos.pt,weight=1)

    h.fill("gen_pos_vxy20_vs_reco_type", vxy=events.GenPos.vxy, reco_type=events.GenPos.matchType, weight=1)
    h.fill("gen_pos_vxy100_vs_reco_type", vxy=events.GenPos.vxy, reco_type=events.GenPos.matchType, weight=1)
    h.fill("gen_pos_pt_vs_reco_type", pt=events.GenPos.pt, reco_type=events.GenPos.matchType, weight=1)

    # ee reco type
    reco_type = ak.to_numpy(np.stack((events.GenEle.matchType,events.GenPos.matchType),axis=1))
    
    mask_reco_R = reco_type == 'R'
    mask_reco_L = reco_type == 'L'
    mask_reco_None = reco_type == 'None'

    mask_reco_RR = ak.values_astype(ak.sum(mask_reco_R, axis=1) == 2, int)
    mask_reco_LL = ak.values_astype(ak.sum(mask_reco_L, axis=1) == 2, int)
    mask_reco_None = ak.values_astype(ak.sum(mask_reco_None, axis=1) > 0, int)
    mask_reco_RL = ak.values_astype((ak.sum(mask_reco_R, axis=1) == 1)&(ak.sum(mask_reco_L, axis=1) == 1), int)

    #print(reco_type)

    vtx_reco_type = mask_reco_LL + 2*mask_reco_RL + 3*mask_reco_RR
    #print(vtx_reco_type)
    
    h.fill("gen_ee_vxy20_vs_reco_type", vxy=events.genEE.vxy, vtype=vtx_reco_type, weight=1)
    h.fill("gen_ee_vxy100_vs_reco_type",vxy=events.genEE.vxy, vtype=vtx_reco_type, weight=1)
    h.fill("gen_ee_pt_vs_reco_type", pt=events.genEE.pt, vtype=vtx_reco_type, weight=1)
    
    #h.fill("gen_vxyCategory_unwgt",vxyCat=events.GenEle.vxyBin,weight=1)
    #h.fill("gen_ptCategory_unwgt",ptCat=events.GenEle.ptBin,weight=1)