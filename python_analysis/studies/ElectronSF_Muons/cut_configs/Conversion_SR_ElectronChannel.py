import numpy as np
import awkward as ak
import sys
sys.path.append("../../../analysisTools/")
import analysisSubroutines as routines
import vector

def cut1(events,info):
    name = "cut1"
    desc = "Pass MET Filters"
    plots = False
    cut = events.METFiltersFailBits == 0
    
    return events[cut], name, desc, plots

def cut2(events,info):
    name = "cut2"
    desc = "HEM Veto"
    plots = False
    if info["year"] == 2018:
        cut = (events.hasHEMjet == 0) & (events.hasHEMelecPF == 0) & (events.hasHEMelecLpt == 0)
        return events[cut], name, desc, plots
    else:
        return events, name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "Single Electron trigger"
    plots = False
    cut = (events.trig.HLT_Ele27_WPTight_Gsf == 1) | (events.trig.HLT_Ele32_WPTight_Gsf == 1) | (events.trig.HLT_Ele35_WPTight_Gsf == 1)

    #print(f'Trigger eff: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "nElectron >= 4"
    plots = False

    nElec = ak.count(events.Electron.pt,axis=1)
    nLptElec = ak.count(events.LptElectron.pt,axis=1)
    cut = (nElec + nLptElec) >= 4

    #print(f'nElectron cut eff: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "n(Prompt Electron) == 2"
    plots = False
    
    nElec_tight = ak.count_nonzero(events.Electron.IDcutMed == 1,axis=1)
    
    cut = nElec_tight >= 2

    #print(f'nTightPromptElectron cut eff: {np.count_nonzero(cut)}/{len(cut)}')
        
    return events[cut], name, desc, plots


def cut6(events,info):
    name = "cut6"
    desc = "Prompt ee pT eta cut"
    plots = False
    
    mask_ele_tight = events.Electron.IDcutMed == 1
    
    #print(events.Electron.pt[mask_ele_tight][ak.count_nonzero(events.Electron.IDcutTight == 1,axis=1)>2]) # pt-sorted already

    # leading and sub-leading
    idx_pt_leading = ak.argmax(ak.argsort(events.Electron.pt[mask_ele_tight],ascending=True),axis=1,keepdims=True)
    idx_pt_subleading = idx_pt_leading + 1
    
    #idx_pt_subleading = ak.argmax(ak.argsort(events.Electron.pt[mask_ele_tight],ascending=True) - ak.argmax(ak.argsort(events.Electron.pt[mask_ele_tight],ascending=True),axis=1), axis=1, keepdims=True)
    #print(idx_pt_subleading)

    sgn_leading = events.Electron.charge[mask_ele_tight][idx_pt_leading]
    sgn_subleading = events.Electron.charge[mask_ele_tight][idx_pt_subleading]

    pt_leading = events.Electron.pt[mask_ele_tight][idx_pt_leading]
    pt_subleading = events.Electron.pt[mask_ele_tight][idx_pt_subleading]

    eta_leading = events.Electron.eta[mask_ele_tight][idx_pt_leading]
    eta_subleading = events.Electron.eta[mask_ele_tight][idx_pt_subleading]

    phi_leading = events.Electron.phi[mask_ele_tight][idx_pt_leading]
    phi_subleading = events.Electron.phi[mask_ele_tight][idx_pt_subleading]

    energy_leading = events.Electron.e[mask_ele_tight][idx_pt_leading]
    energy_subleading = events.Electron.e[mask_ele_tight][idx_pt_subleading]

#    dict_prompt_vtx = {'e1': {'pt': ak.flatten(pt_leading), 'eta': ak.flatten(eta_leading), 'phi': ak.flatten(phi_leading)},\
#                       'e2': {'pt': ak.flatten(pt_subleading), 'eta': ak.flatten(eta_subleading), 'phi': ak.flatten(phi_subleading)}}
#    events.__setitem__("prompt_vtx",dict_prompt_vtx)

    events.__setitem__("prompt_e1_pt",ak.flatten(pt_leading))
    events.__setitem__("prompt_e1_eta",ak.flatten(eta_leading))
    events.__setitem__("prompt_e1_phi",ak.flatten(phi_leading))
    events.__setitem__("prompt_e1_e",ak.flatten(energy_leading))

    events.__setitem__("prompt_e2_pt",ak.flatten(pt_subleading))
    events.__setitem__("prompt_e2_eta",ak.flatten(eta_subleading))
    events.__setitem__("prompt_e2_phi",ak.flatten(phi_subleading))
    events.__setitem__("prompt_e2_e",ak.flatten(energy_subleading))

    # leading pT cut with trigger plateau value
    trig_ele27 = (events.trig.HLT_Ele27_WPTight_Gsf == 1) & (events.trig.HLT_Ele32_WPTight_Gsf == 0) & (events.trig.HLT_Ele35_WPTight_Gsf == 0)
    trig_ele32 = (events.trig.HLT_Ele32_WPTight_Gsf == 1) & (events.trig.HLT_Ele35_WPTight_Gsf == 0)
    trig_ele35 = (events.trig.HLT_Ele35_WPTight_Gsf == 1)

    cut_pt_leading27 = (pt_leading > 30) & trig_ele27
    cut_pt_leading32 = (pt_leading > 35) & trig_ele32
    cut_pt_leading35 = (pt_leading > 38) & trig_ele35
    
#    cut_pt_leading27 = (pt_leading > 37) & trig_ele27
#    cut_pt_leading32 = (pt_leading > 42) & trig_ele32
#    cut_pt_leading35 = (pt_leading > 45) & trig_ele35
    
    cut_pt_leading = cut_pt_leading27 | cut_pt_leading32 | cut_pt_leading35 
    cut_pt_subleading = (pt_subleading > 15)
#    cut_pt_subleading = (pt_subleading > 10)

    cut_pt = ak.flatten(cut_pt_leading) & ak.flatten(cut_pt_subleading)

    # eta cut
    cut_eta = ak.flatten(np.abs(eta_leading) < 2.4) & ak.flatten(np.abs(eta_subleading) < 2.4)

    # OS
    cut_OS = ak.flatten(sgn_leading + sgn_subleading) == 0
    
    cut = cut_pt & cut_eta & cut_OS

    #print(f'prompt ee pT & eta & OS cut eff: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots

def cut7(events,info):
    name = "cut7"
    desc = "sel_vtx == low pT collections"
    plots = False

    if len(events) == 0:
        cut = []
    else:
        cut = (events.sel_vtx.e1_typ == "L") & (events.sel_vtx.e2_typ == "L")
    #    cut = (events.sel_vtx.e1_typ == "R") & (events.sel_vtx.e2_typ == "R")
        
        #print(f'selvtx == LL cut eff: {np.count_nonzero(cut)}/{len(cut)}')
    
        # Save mee and meeee 
        #mee = routines.getInvariantMass([events.prompt_e1_pt, events.prompt_e2_pt], [events.prompt_e1_eta, events.prompt_e2_eta], [events.prompt_e1_phi, events.prompt_e2_phi])
        #events.__setitem__("prompt_mee",mee)
    
        behavior = vector.backends.awkward.behavior
        vec_prompt_l1 = ak.zip(
            {
                "pt":events.prompt_e1_pt,
                "eta":events.prompt_e1_eta,
                "phi":events.prompt_e1_phi,
                "energy":events.prompt_e1_e
            },
            with_name="Momentum4D",
            behavior=behavior
        )
    
        vec_prompt_l2 = ak.zip(
            {
                "pt":events.prompt_e2_pt,
                "eta":events.prompt_e2_eta,
                "phi":events.prompt_e2_phi,
                "energy":events.prompt_e2_e
            },
            with_name="Momentum4D",
            behavior=behavior
        )
    
        vec_displaced_ee = ak.zip(
            {
                "pt":events.sel_vtx.pt,
                "eta":events.sel_vtx.eta,
                "phi":events.sel_vtx.phi,
                "mass":events.sel_vtx.m
            },
            with_name="Momentum4D",
            behavior=behavior
        )
    
        vec_displaced_ee_refit = ak.zip(
            {
                "pt":events.sel_vtx.refit_pt,
                "eta":events.sel_vtx.refit_eta,
                "phi":events.sel_vtx.refit_phi,
                "mass":events.sel_vtx.refit_m
            },
            with_name="Momentum4D",
            behavior=behavior
        )
        #behavior.update(ak._util.copy_behaviors("Momentum4D", "LorentzVector", behavior))
    
        '''
        vec_prompt_e1 = vector.obj(pt=events.prompt_e1_pt, phi=events.prompt_e1_phi, eta=events.prompt_e1_eta, E=events.prompt_e1_e)
    
        print(vec_prompt_e1.eta)
        
        vec_prompt_e2 = vector.obj(pt=events.prompt_e2_pt, phi=events.prompt_e2_phi, eta=events.prompt_e2_eta, E=events.prompt_e2_e)
        '''
        
        vec_prompt_ee = vec_prompt_l1 + vec_prompt_l2
        vec_eeee = vec_prompt_l1 + vec_prompt_l2 + vec_displaced_ee
        vec_eeee_refit = vec_prompt_l1 + vec_prompt_l2 + vec_displaced_ee_refit
    
        '''
        print('e1 eta', vec_prompt_e1.eta)
        print('e2 eta', vec_prompt_e2.eta)
        print('ee eta', vec_prompt_ee.eta)
    
        print('ee mass', vec_prompt_ee.mass)
        print('mee', mee)
        '''
        #vec_e1 = vector.obj(pt=events.prompt_e1_pt, eta=events.prompt_e1_eta, phi=events.prompt_e1_phi)
        #vec_e2 = vector.obj(pt=events.prompt_e2_pt, eta=events.prompt_e2_eta, phi=events.prompt_e2_phi)
    
        events.prompt_e1_px = events.prompt_e1_pt * np.cos(events.prompt_e1_phi)
        events.prompt_e1_py = events.prompt_e1_pt * np.sin(events.prompt_e1_phi)
        events.prompt_e1_pz = events.prompt_e1_pt * np.sinh(events.prompt_e1_eta)
    
        events.prompt_e2_px = events.prompt_e2_pt * np.cos(events.prompt_e2_phi)
        events.prompt_e2_py = events.prompt_e2_pt * np.sin(events.prompt_e2_phi)
        events.prompt_e2_pz = events.prompt_e2_pt * np.sinh(events.prompt_e2_eta)
    
        events.prompt_ee_px = events.prompt_e1_px + events.prompt_e2_px
        events.prompt_ee_py = events.prompt_e1_py + events.prompt_e2_py
        events.prompt_ee_pz = events.prompt_e1_pz + events.prompt_e2_pz
    
        ee_prompt_pt = np.sqrt(events.prompt_ee_px**2 + events.prompt_ee_py**2)
        ee_prompt_eta = np.arcsinh(events.prompt_ee_pz/ee_prompt_pt)
        ee_prompt_phi = np.arccos(events.prompt_ee_px/ee_prompt_pt)
    
        #print('manual ee eta', ee_prompt_eta)
    
        #ee_prompt_pt = events.prompt_e1_pt + events.prompt_e2_pt
        #ee_prompt_eta = events.prompt_e1_eta + events.prompt_e2_eta
        #ee_prompt_phi = events.prompt_e1_phi + events.prompt_e2_phi
    
        ee_disp_pt = events.sel_vtx.pt
        ee_disp_eta = events.sel_vtx.eta
        ee_disp_phi = events.sel_vtx.phi
    
        ee_disp_pt_refit = events.sel_vtx.refit_pt
        ee_disp_eta_refit = events.sel_vtx.refit_eta
        ee_disp_phi_refit = events.sel_vtx.refit_phi
        
        #meeee = routines.getInvariantMass([ee_prompt_pt,ee_disp_pt], [ee_prompt_eta,ee_disp_eta], [ee_prompt_phi,ee_disp_phi])
        #meeee_refit = routines.getInvariantMass([ee_prompt_pt,ee_disp_pt_refit], [ee_prompt_eta,ee_disp_eta_refit], [ee_prompt_phi,ee_disp_phi_refit])

        pt_prompt_lead = vec_prompt_l1.pt
        events.__setitem__("pt_prompt_lead",pt_prompt_lead)
        pt_prompt_sublead = vec_prompt_l2.pt
        events.__setitem__("pt_prompt_sublead",pt_prompt_sublead)
        
        pt_prompt_2l = vec_prompt_ee.pt
        events.__setitem__("pt_prompt_2l",pt_prompt_2l)
        
        mee = vec_prompt_ee.mass
        events.__setitem__("prompt_m2l",mee)
    
        meeee = vec_eeee.mass
        events.__setitem__("m4l",meeee)
    
        meeee_refit = vec_eeee_refit.mass
        events.__setitem__("m4l_refit",meeee_refit)
        
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "dR(prompt, disp) > 0.05"
    plots = True

    
    if len(events) == 0:
        cut = []
    else:
        cut = routines.deltaR(events.sel_vtx.e1.eta, events.sel_vtx.e1.phi, events.prompt_e1_eta, events.prompt_e1_phi) > 0.05
        cut = cut & (routines.deltaR(events.sel_vtx.e1.eta, events.sel_vtx.e1.phi, events.prompt_e2_eta, events.prompt_e2_phi) > 0.05)
        cut = cut & (routines.deltaR(events.sel_vtx.e2.eta, events.sel_vtx.e2.phi, events.prompt_e1_eta, events.prompt_e1_phi) > 0.05)
        cut = cut & (routines.deltaR(events.sel_vtx.e2.eta, events.sel_vtx.e2.phi, events.prompt_e2_eta, events.prompt_e2_phi) > 0.05)
        
        #print(f'overlap removal btw prompt & displaced: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots


def cut9(events,info):
    name = "cut9"
    desc = "M(ee) refit < 0.1"
    plots = True

#    nMuon_passID = ak.count_nonzero(events.Muon.IDcutMediumPrompt == 1,axis=1)
    #nMuon_passID = ak.count_nonzero(events.Muon.IDcutTight == 1,axis=1)
    
    cut = events.sel_vtx.refit_m < 0.1
    
    return events[cut], name, desc, plots


def cut10(events,info):
    name = "cut10"
    desc = "80 < M(eeee) < 100"
    plots = True
    
    if len(events) == 0:
        cut = []
    else:
        cut = (events.m4l_refit > 80) & (events.m4l_refit < 100)
    
    print(f'M(eeee) cut: {np.count_nonzero(cut)}/{len(cut)}')
        
    return events[cut], name, desc, plots

def cut11(events,info):
    name = "cut11"
    desc = "50 < M(ee) prompt"
    plots = True

    
    if len(events) == 0:
        cut = []
    else:
        cut = (events.prompt_m2l > 50)
    
    #print(f'M(ee) prompt cut: {np.count_nonzero(cut)}/{len(cut)}')

    return events[cut], name, desc, plots

def cut12(events,info):
    name = "cut12"
    desc = "M(ee) prompt < 80"
    plots = True

    
    if len(events) == 0:
        cut = []
    else:
        cut = (events.prompt_m2l < 80)
    
    #print(f'M(ee) prompt cut: {np.count_nonzero(cut)}/{len(cut)}')

    return events[cut], name, desc, plots
