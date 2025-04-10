import numpy as np
import awkward as ak
from histoBinning import myHisto

def make_histograms(info):
    h = myHisto()

    # electron
    h.make("sel_e1_dxy",'dxy')
    h.make("sel_e1_dxy_refit",'dxy')
    h.make("sel_e1_vz",'vz')

    h.make("sel_e2_dxy",'dxy')
    h.make("sel_e2_dxy_refit",'dxy')
    h.make("sel_e2_vz",'vz')

    h.make("sel_e1e2_dxy_refit", 'dxy')

    h.make("PVx", 'pvx')
    h.make("PVy", 'pvy')


    ## Refitting of vtx track ##
    h.make('sel_vtx_mass_refit','vtx_mass_refit')
    h.make('sel_vtx_mass_low_refit','mass_low_refit')
    h.make('sel_vtx_pt_refit','vtx_pt_refit')
    h.make('sel_vtx_eta_refit','eta_refit')
    h.make('sel_vtx_phi_refit','phi_refit')
    h.make('sel_vtx_dR_refit','dR_refit')

    h.make('sel_vtx_vxy1','vxy1')
    h.make('sel_vtx_vxy10','vxy10')
    h.make('sel_vtx_vxy100','vxy100')
    
    h.make('sel_vtx_vxy1_fromPV','vxy1')
    h.make('sel_vtx_vxy10_fromPV','vxy10')
    h.make('sel_vtx_vxy100_fromPV','vxy100')

    h.make('sel_vtx_vxy_fromPV','vxy_SFbin')
    h.make('sel_vtx_vxy_fromPV_coarse','vxy_SFbin_coarse')
    h.make('sel_vtx_vxy_fromPV_to5','vxy_SFbin_to5')

    h.make('sel_vtx_vxy_fromPV_coarse_vs_mass_low_refit','vxy_SFbin_coarse','mass_low_refit')
    h.make('sel_vtx_vxy_fromPV_to5_vs_mass_low_refit','vxy_SFbin_to5','mass_low_refit')
    h.make('sel_vtx_vxy_fromPV_vs_mass_low_refit','vxy_SFbin','mass_low_refit')

    h.make('prompt_vtx_mass','vtx_mass_high')
    h.make('mass_fourlepton','vtx_mass_high')

    h.make('mass_fourlepton_refit','vtx_mass_high')

    h.make("sel_e1_log10dxydz","log10dEtadPhi")
    h.make("sel_e2_log10dxydz","log10dEtadPhi")
    h.make("min_sel_log10dxydz","log10dEtadPhi")

    h.make('sel_vtx_vx_vy','vx','vy')
    
    h.make("pt_prompt_lead","vtx_pt")
    h.make("pt_prompt_sublead","vtx_pt")
    h.make("pt_prompt_2l","vtx_pt")
    h.make('sel_vtx_pt_vs_prompt_pt','vtx_pt_refit','vtx_pt')
    h.make('sel_vtx_vxy_vs_prompt_pt','vxy_SFbin','vtx_pt')
    
    h.make('sel_vtx_vxy_vs_pt','vxy_SFbin','vtx_pt_refit')
    h.make('sel_vtx_vxy_vs_dR','vxy_SFbin','dR_refit')

    h.make('nPV','nPV')
    
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
    h.fill("sel_e1_dxy",dxy=sel_vtx.e1.dxy,weight=wgt)
    h.fill("sel_e1_dxy_refit", dxy=np.abs(sel_vtx.e1.refit_dxy), weight=wgt)
    h.fill("sel_e1_vz",vz=sel_vtx.e1.vz,weight=wgt)
    
    h.fill("sel_e2_dxy",dxy=sel_vtx.e2.dxy,weight=wgt)
    h.fill("sel_e2_dxy_refit", dxy=np.abs(sel_vtx.e2.refit_dxy), weight=wgt)
    h.fill("sel_e2_vz",vz=sel_vtx.e2.vz,weight=wgt)

    if info["type"] == "signal" or info["type"] == "bkg":
        h.fill("sel_e1e2_dxy_refit", dxy=np.concatenate([np.abs(sel_vtx.e1.refit_dxy), np.abs(sel_vtx.e2.refit_dxy)]), weight=np.concatenate([wgt, wgt]))
    else:
        h.fill("sel_e1e2_dxy_refit", dxy=np.concatenate([np.abs(sel_vtx.e1.refit_dxy), np.abs(sel_vtx.e2.refit_dxy)]), weight=wgt)

    h.fill("PVx",pvx=events.PV.x,weight=wgt)
    h.fill("PVy",pvy=events.PV.y,weight=wgt)
    
    ## Refit
    h.fill('sel_vtx_mass_refit',mass_refit=sel_vtx.refit_m,weight=wgt)
    h.fill('sel_vtx_mass_low_refit',mass_low_refit=sel_vtx.refit_m,weight=wgt)
    h.fill('sel_vtx_pt_refit',pt_refit=sel_vtx.refit_pt,weight=wgt)
    h.fill('sel_vtx_eta_refit',eta_refit=sel_vtx.refit_eta,weight=wgt)
    h.fill('sel_vtx_phi_refit',phi_refit=sel_vtx.refit_phi,weight=wgt)
    h.fill('sel_vtx_dR_refit',dR_refit=sel_vtx.refit_dR,weight=wgt)

    h.fill('sel_vtx_vxy1',vxy=np.sqrt((sel_vtx.vx)**2+(sel_vtx.vy)**2),weight=wgt)
    h.fill('sel_vtx_vxy10',vxy=np.sqrt((sel_vtx.vx)**2+(sel_vtx.vy)**2),weight=wgt)
    h.fill('sel_vtx_vxy100',vxy=np.sqrt((sel_vtx.vx)**2+(sel_vtx.vy)**2),weight=wgt)

    h.fill('sel_vtx_vxy1_fromPV',vxy=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),weight=wgt)
    h.fill('sel_vtx_vxy10_fromPV',vxy=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),weight=wgt)
    h.fill('sel_vtx_vxy100_fromPV',vxy=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),weight=wgt)

    h.fill('sel_vtx_vxy_fromPV',vxy_SFbin=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),weight=wgt)
    h.fill('sel_vtx_vxy_fromPV_coarse',vxy_SFbin_coarse=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),weight=wgt)

    h.fill('sel_vtx_vxy_fromPV_to5',vxy_SFbin_coarse=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),weight=wgt)
    
    h.fill('prompt_vtx_mass',mass=events.prompt_m2l,weight=wgt)
    h.fill('mass_fourlepton',mass=events.m4l,weight=wgt)

    h.fill('mass_fourlepton_refit',mass=events.m4l_refit,weight=wgt)

    h.fill('sel_vtx_vx_vy',vx=sel_vtx.vx,vy=sel_vtx.vy,weight=wgt)

    h.fill('nPV', nPV=events.numPV,weight=wgt)

    h.fill('sel_vtx_vxy_fromPV_coarse_vs_mass_low_refit',vxy_SFbin_coarse=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),mass_low_refit=sel_vtx.refit_m,weight=wgt)
    h.fill('sel_vtx_vxy_fromPV_to5_vs_mass_low_refit',vxy_SFbin_coarse=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),mass_low_refit=sel_vtx.refit_m,weight=wgt)
    h.fill('sel_vtx_vxy_fromPV_vs_mass_low_refit',vxy_SFbin=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),mass_low_refit=sel_vtx.refit_m,weight=wgt)
    

    h.fill("pt_prompt_lead",pt=events.pt_prompt_lead,weight=wgt)
    h.fill("pt_prompt_sublead",pt=events.pt_prompt_sublead,weight=wgt)
    h.fill("pt_prompt_2l",pt=events.pt_prompt_2l,weight=wgt)
    h.fill('sel_vtx_pt_vs_prompt_pt',pt_refit=sel_vtx.refit_pt,pt=events.pt_prompt_2l,weight=wgt)
    h.fill('sel_vtx_vxy_vs_prompt_pt',vxy_SFbin=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),pt=events.pt_prompt_2l,weight=wgt)
    
    h.fill('sel_vtx_vxy_vs_pt',vxy_SFbin=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),pt_refit=sel_vtx.refit_pt,weight=wgt)
    h.fill('sel_vtx_vxy_vs_dR',vxy_SFbin=np.sqrt((sel_vtx.vx-events.PV.x)**2+(sel_vtx.vy-events.PV.y)**2),dR_refit=sel_vtx.refit_dR,weight=wgt)

    h.fill("sel_e1_log10dxydz",log10dEtadPhi=np.log10(np.abs(sel_vtx.e1.dxy/sel_vtx.e1.dz)),weight=wgt)
    h.fill("sel_e2_log10dxydz",log10dEtadPhi=np.log10(np.abs(sel_vtx.e2.dxy/sel_vtx.e2.dz)),weight=wgt)
    h.fill("min_sel_log10dxydz",log10dEtadPhi=np.minimum(np.log10(np.abs(sel_vtx.e1.dxy/sel_vtx.e1.dz)), 
                                                         np.log10(np.abs(sel_vtx.e2.dxy/sel_vtx.e2.dz))),weight=wgt)
