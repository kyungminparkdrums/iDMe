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
    desc = "Single Muon trigger"
    plots = False
    cut = (events.trig.HLT_IsoMu27 == 1)

    print(f'Trigger eff: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "nLptElectron >= 2"
    plots = False

    nLptElec = ak.count(events.LptElectron.pt,axis=1)
    cut = nLptElec >= 2

    print(f'nLptElec cut eff: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "nMuon == 2"
    plots = False

    nMuon = ak.count(events.Muon.pt,axis=1)
    cut = nMuon >= 2

    print(f'nMuon cut eff: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots

def cut6(events,info):
    name = "cut6"
    desc = "nMuon >= 2 with PASS ID"
    plots = False

#    nMuon_passID = ak.count_nonzero(events.Muon.IDcutMediumPrompt == 1,axis=1)
    #nMuon_passID = ak.count_nonzero(events.Muon.IDcutTight == 1,axis=1)
    nMuon_passID = ak.count_nonzero(events.Muon.IDcutMedium == 1,axis=1)

    cut = nMuon_passID >= 2

    print(f'nMuon cut eff: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots

def cut7(events,info):
    name = "cut7"
    desc = "Prompt di-muon pT eta cut"
    plots = False
    
    #mask_mu_tight = events.Muon.IDcutTight == 1
    mask_mu_tight = events.Muon.IDcutMedium == 1
    
    # leading and sub-leading
    idx_pt_leading = ak.argmax(ak.argsort(events.Muon.pt[mask_mu_tight],ascending=True),axis=1,keepdims=True)
    idx_pt_subleading = idx_pt_leading + 1

    sgn_leading = events.Muon.charge[mask_mu_tight][idx_pt_leading]
    sgn_subleading = events.Muon.charge[mask_mu_tight][idx_pt_subleading]

    pt_leading = events.Muon.pt[mask_mu_tight][idx_pt_leading]
    pt_subleading = events.Muon.pt[mask_mu_tight][idx_pt_subleading]

    eta_leading = events.Muon.eta[mask_mu_tight][idx_pt_leading]
    eta_subleading = events.Muon.eta[mask_mu_tight][idx_pt_subleading]

    phi_leading = events.Muon.phi[mask_mu_tight][idx_pt_leading]
    phi_subleading = events.Muon.phi[mask_mu_tight][idx_pt_subleading]

    energy_leading = events.Muon.e[mask_mu_tight][idx_pt_leading]
    energy_subleading = events.Muon.e[mask_mu_tight][idx_pt_subleading]

    # set event quantity
    events.__setitem__("prompt_mu1_pt",ak.flatten(pt_leading))
    events.__setitem__("prompt_mu1_eta",ak.flatten(eta_leading))
    events.__setitem__("prompt_mu1_phi",ak.flatten(phi_leading))
    events.__setitem__("prompt_mu1_e",ak.flatten(energy_leading))

    events.__setitem__("prompt_mu2_pt",ak.flatten(pt_subleading))
    events.__setitem__("prompt_mu2_eta",ak.flatten(eta_subleading))
    events.__setitem__("prompt_mu2_phi",ak.flatten(phi_subleading))
    events.__setitem__("prompt_mu2_e",ak.flatten(energy_subleading))

    # leading pT cut with trigger plateau value
    cut_pt_leading = pt_leading > 28
    #cut_pt_leading = pt_leading > 30
    cut_pt_subleading = (pt_subleading > 15)

    cut_pt = ak.flatten(cut_pt_leading) & ak.flatten(cut_pt_subleading)

    # eta cut
    cut_eta = ak.flatten(np.abs(eta_leading) < 2.4) & ak.flatten(np.abs(eta_subleading) < 2.4)

    # OS
    cut_OS = ak.flatten(sgn_leading + sgn_subleading) == 0
    
    cut = cut_pt & cut_eta & cut_OS

    print(f'prompt di-muon pT & eta & OS cut eff: {np.count_nonzero(cut)}/{len(cut)}')
    
    return events[cut], name, desc, plots

def cut8(events,info):
    name = "cut8"
    desc = "sel_vtx == low pT collections"
    plots = False

    if len(events) == 0:
        cut = []
    else:
        cut = (events.sel_vtx.e1_typ == "L") & (events.sel_vtx.e2_typ == "L")
        
        print(f'selvtx == LL cut eff: {np.count_nonzero(cut)}/{len(cut)}')
    
        # Save mee and meeee
        behavior = vector.backends.awkward.behavior
        vec_prompt_l1 = ak.zip(
            {
                "pt":events.prompt_mu1_pt,
                "eta":events.prompt_mu1_eta,
                "phi":events.prompt_mu1_phi,
                "energy":events.prompt_mu1_e
            },
            with_name="Momentum4D",
            behavior=behavior
        )
    
        vec_prompt_l2 = ak.zip(
            {
                "pt":events.prompt_mu2_pt,
                "eta":events.prompt_mu2_eta,
                "phi":events.prompt_mu2_phi,
                "energy":events.prompt_mu2_e
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
        
        vec_prompt_mumu = vec_prompt_l1 + vec_prompt_l2
        vec_4l = vec_prompt_l1 + vec_prompt_l2 + vec_displaced_ee
        vec_4l_refit = vec_prompt_l1 + vec_prompt_l2 + vec_displaced_ee_refit
    
        events.prompt_mu1_px = events.prompt_mu1_pt * np.cos(events.prompt_mu1_phi)
        events.prompt_mu1_py = events.prompt_mu1_pt * np.sin(events.prompt_mu1_phi)
        events.prompt_mu1_pz = events.prompt_mu1_pt * np.sinh(events.prompt_mu1_eta)
    
        events.prompt_mu2_px = events.prompt_mu2_pt * np.cos(events.prompt_mu2_phi)
        events.prompt_mu2_py = events.prompt_mu2_pt * np.sin(events.prompt_mu2_phi)
        events.prompt_mu2_pz = events.prompt_mu2_pt * np.sinh(events.prompt_mu2_eta)
    
        events.prompt_mumu_px = events.prompt_mu1_px + events.prompt_mu2_px
        events.prompt_mumu_py = events.prompt_mu1_py + events.prompt_mu2_py
        events.prompt_mumu_pz = events.prompt_mu1_pz + events.prompt_mu2_pz
    
        mumu_prompt_pt = np.sqrt(events.prompt_mumu_px**2 + events.prompt_mumu_py**2)
        mumu_prompt_eta = np.arcsinh(events.prompt_mumu_pz/mumu_prompt_pt)
        mumu_prompt_phi = np.arccos(events.prompt_mumu_px/mumu_prompt_pt)

        ee_disp_pt = events.sel_vtx.pt
        ee_disp_eta = events.sel_vtx.eta
        ee_disp_phi = events.sel_vtx.phi
    
        ee_disp_pt_refit = events.sel_vtx.refit_pt
        ee_disp_eta_refit = events.sel_vtx.refit_eta
        ee_disp_phi_refit = events.sel_vtx.refit_phi

        pt_prompt_lead = vec_prompt_l1.pt
        events.__setitem__("pt_prompt_lead",pt_prompt_lead)
        pt_prompt_sublead = vec_prompt_l2.pt
        events.__setitem__("pt_prompt_sublead",pt_prompt_sublead)
        
        pt_prompt_2l = vec_prompt_mumu.pt
        events.__setitem__("pt_prompt_2l",pt_prompt_2l)
        
        m2l = vec_prompt_mumu.mass
        events.__setitem__("prompt_m2l",m2l)
    
        m4l = vec_4l.mass
        events.__setitem__("m4l",m4l)
    
        m4l_refit = vec_4l_refit.mass
        events.__setitem__("m4l_refit",m4l_refit)
        
    return events[cut], name, desc, plots


def cut9(events,info):
    name = "cut9"
    desc = "M(ee) refit < 0.1"
    plots = True
    
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
        #cut = (events.prompt_mee < 80)
    
    #print(f'M(ee) prompt cut: {np.count_nonzero(cut)}/{len(cut)}')

    return events[cut], name, desc, plots

def cut12(events,info):
    name = "cut12"
    desc = "M(ee) prompt < 80"
    plots = True

    
    if len(events) == 0:
        cut = []
    else:
        #cut = (events.prompt_m2l < 90)
        cut = (events.prompt_m2l < 80)
    
    #print(f'M(ee) prompt cut: {np.count_nonzero(cut)}/{len(cut)}')

    return events[cut], name, desc, plots
