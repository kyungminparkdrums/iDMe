import correctionlib
import numpy as np
import awkward as ak

from coffea.lookup_tools import extractor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory



##########
## PU SF
###########
def get_sf_PU(IOV, jsonPath, nTrueInt, type):
    # Year
    if (IOV=='2018'): 
        key = 'Collisions18_UltraLegacy_goldenJSON'
    elif (IOV=='2017'): 
        key = 'Collisions17_UltraLegacy_goldenJSON'
    elif (IOV=='2016'):  # postVFP
        key = 'Collisions16_UltraLegacy_goldenJSON'
    elif (IOV=='2016APV'):  # preVFP
        key = 'Collisions16_UltraLegacy_goldenJSON'
    
    evaluator = correctionlib.CorrectionSet.from_file(f'{jsonPath}/puWeights.json')

    # Nominal/Up/Down
    sf = evaluator[key].evaluate(nTrueInt, type) # type = 'nominal'/'up'/'down'

    return sf

##########
## Electron SF
###########

def get_sf_elePFid(IOV, jsonPath, pt, eta, type):
    ## Reference:
    ##   - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
    ##   - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaSFJSON
    ## json files from: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM

    # Year
    year = {
        "2016APV" : "2016preVFP",
        "2016"    : "2016postVFP",
        "2017"    : "2017",
        "2018"    : "2018",
    }
    num = ak.num(pt)
    evaluator = correctionlib.CorrectionSet.from_file(f'{jsonPath}/electron.json')

    # Type
    if type == 'nominal':
        corrType = 'sf'
    elif type == 'up':
        corrType = 'sfup'
    elif type == 'down':
        corrType = 'sfdown'
    
    ## EGM POG doesn't provide SF under 10 GeV; set it to 1 here and apply later
    mask = pt > 10
    pt = ak.where(mask, pt, 12)
    
    sf = evaluator["UL-Electron-ID-SF"].evaluate(year[IOV], corrType, "Loose",
                                                 np.array(ak.flatten(eta)),
                                                 np.array(ak.flatten(pt)))
    sf = ak.where(np.array(ak.flatten(~mask)), 1, sf)

    return ak.unflatten(sf, ak.num(pt))

def get_sf_elePFreco(IOV, jsonPath, pt, eta, type):
    ## Reference:
    ##   - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
    ##   - https://twiki.cern.ch/twiki/bin/view/CMS/EgammaSFJSON
    ## json files from: https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/EGM

    # Year
    year = {
        "2016APV" : "2016preVFP",
        "2016"    : "2016postVFP",
        "2017"    : "2017",
        "2018"    : "2018",
    }
    num = ak.num(pt)
    evaluator = correctionlib.CorrectionSet.from_file(f'{jsonPath}/electron.json')

    # Type
    if type == 'nominal':
        corrType = 'sf'
    elif type == 'up':
        corrType = 'sfup'
    elif type == 'down':
        corrType = 'sfdown'

    # Reco [10,20]
    mask_below = (pt < 20) & (pt > 10)
    pt_below = ak.where(mask_below, pt, 15)
    sf_below = evaluator["UL-Electron-ID-SF"].evaluate(year[IOV], corrType, "RecoBelow20",
                                                 np.array(ak.flatten(eta)),
                                                 np.array(ak.flatten(pt_below)))
    sf_below = ak.where(np.array(ak.flatten(~mask_below)), 0, sf_below)

    # Reco [20,inf]
    mask_above = pt >= 20
    pt_above = ak.where(mask_above, pt, 25)
    sf_above = evaluator["UL-Electron-ID-SF"].evaluate(year[IOV], corrType, "RecoAbove20",
                                                 np.array(ak.flatten(eta)),
                                                 np.array(ak.flatten(pt_above)))
    sf_above = ak.where(np.array(ak.flatten(~mask_above)), 0, sf_above)

    # SF
    sf = sf_above + sf_below
    mask_zeroSF = (sf == 0) # for [-inf, 10] GeV bins

    sf = ak.where(np.array(mask_zeroSF), 1, sf)
    
    return ak.unflatten(sf, ak.num(pt))

def get_sf_elePF5to10(IOV, jsonPath, pt, eta, type):
    evaluator = correctionlib.CorrectionSet.from_file(f'{jsonPath}/electronPF5to10_ptbin.json')
    
    mask = (pt > 5.) & (pt < 10.)
    pt = ak.where(mask, pt, 7.5)

    sf = evaluator['PFelectron5to10GeV'].evaluate(np.array(ak.flatten(pt)), type) # type = 'nominal'/'up'/'down'
    sf = ak.where(np.array(ak.flatten(~mask)), 0, sf)

    return ak.unflatten(sf, ak.num(pt))

def get_sf_elePF(IOV, jsonPath, pt, eta, type):
    # from [10, inf] bins as provided by EGM POG
    sf_elePFreco = get_sf_elePFreco(IOV, jsonPath, pt, eta, type)
    sf_elePFid = get_sf_elePFid(IOV, jsonPath, pt, eta, type)

    mask_ptAbove10 = pt >= 10

    sf_above10 = (sf_elePFreco * sf_elePFid)
    sf_above10 = ak.where(np.array(ak.flatten(~mask_ptAbove10)), 0, ak.flatten(sf_above10))
    
    # for [5,10] GeV, put 1 for now
    sf_below10 = get_sf_elePF5to10(IOV, jsonPath, pt, eta, type)
    sf_below10 = ak.where(np.array(ak.flatten(mask_ptAbove10)), 0, ak.flatten(sf_below10))

    # PF electron SF across all pT spectrum
    sf = sf_above10 + sf_below10
    sf = ak.unflatten(sf, ak.num(pt))
    
    return sf

def get_sf_eleLpt(jsonPath, dxy, type):
    evaluator = correctionlib.CorrectionSet.from_file(f'{jsonPath}/lptElectron.json')
    sf = evaluator['LptElectron'].evaluate(np.array(ak.flatten(dxy)), type) # type = 'nominal'/'up'/'down'

    print(sf)
    
    return ak.unflatten(sf, ak.num(dxy))

def get_sf_ee_vtx(IOV, jsonPath, vxy, vtx_type, type):
    type_list = ['LL', 'LR', 'RR']

    evaluator = {}
    for typ in type_list:
        evaluator[typ] = correctionlib.CorrectionSet.from_file(f'{jsonPath}/ee_vtx_{typ}_{IOV}.json')
        
    mask_vxy = (vxy > 20)
    vxy = ak.where(mask_vxy, 15, vxy)

    mask_type = {}
    sf_type = {}
    for typ in type_list:
        mask_type[typ] = vtx_type == typ
        sf_type[typ] = evaluator[typ][f'ee_vtx_{typ}'].evaluate(np.array(ak.flatten(vxy)), type) # type = 'nominal'/'up'/'down'
        sf_type[typ] = ak.where(ak.flatten(~mask_type[typ]), 0, sf_type[typ])
    sf = sf_type['LL'] + sf_type['LR'] + sf_type['RR']

    sf = ak.unflatten(sf, ak.num(vxy))
    
    return sf

##########
## Trigger
##########
def get_trigger_MC(IOV, jsonPath, MET, isMC, type):
    if isMC:
        evaluator = correctionlib.CorrectionSet.from_file(f'{jsonPath}/trig_MC_{IOV}.json')
    else:
        evaluator = correctionlib.CorrectionSet.from_file(f'{jsonPath}/trig_data_{IOV}.json')
    
    mask_MET = (MET > 1000)
    MET = ak.where(mask_MET, 800, MET)

    sf = evaluator[typ][f'trigger_MC'].evaluate(np.array(ak.flatten(MET)), type) # type = 'nominal'/'up'/'down' # trigger_MC naming convention bug; this name is also applied to data

    print(sf)
    
    return sf



##########
## b-tagging SF
###########
def get_btagSF(jsonPath, pt, eta, truth, type):
    evaluator = correctionlib.CorrectionSet.from_file(f'{jsonPath}/btagging.json')

    if type == 'nominal':
        bc_corrType = 'central'
        light_corrType = 'central'
    elif type == 'bc_correlated_up':
        bc_corrType = 'up_correlated'
        light_corrType = 'central'
    elif type == 'light_correlated_up':
        bc_corrType = 'central'
        light_corrType = 'up_correlated'
    elif type == 'bc_uncorrelated_up':
        bc_corrType = 'up_uncorrelated'
        light_corrType = 'central'
    elif type == 'light_uncorrelated_up':
        bc_corrType = 'central'
        light_corrType = 'up_uncorrelated'
    elif type == 'bc_correlated_down':
        bc_corrType = 'down_correlated'
        light_corrType = 'central'
    elif type == 'light_correlated_down':
        bc_corrType = 'central'
        light_corrType = 'down_correlated'
    elif type == 'bc_uncorrelated_down':
        bc_corrType = 'down_uncorrelated'
        light_corrType = 'central'
    elif type == 'light_uncorrelated_down':
        bc_corrType = 'central'
        light_corrType = 'down_uncorrelated'
    
    abseta = np.array(ak.flatten(np.abs(eta)))
    abseta_mask = abseta > 2.4
    abseta[abseta_mask] = 2.4

    #print(abseta_mask)

    #print(abseta_mask)
    truth = np.array(ak.flatten(truth))

    l_mask = (truth == 0)
    bc_flav =  np.array(ak.where(l_mask, 4, truth))
    l_flav =  np.array(ak.zeros_like(truth))
    
    bc_sf = evaluator["deepJet_comb"].evaluate(bc_corrType, 'M', bc_flav, abseta, np.array(ak.flatten(pt)))
    light_sf = evaluator["deepJet_incl"].evaluate(light_corrType, 'M', l_flav, abseta, np.array(ak.flatten(pt)))
    
    sf = ak.where(l_mask, light_sf, bc_sf)
    sf = ak.where(abseta_mask, 1, sf)
    #print(sf)
    
    return ak.unflatten(sf, ak.num(pt))

def apply_btagSF(iov, jsonPath, events, type):
    mask_isTagged = events.PFJet.passMedID # the ones that don't MED ID will be used for SF derivation effectively

    truth = np.array(ak.flatten(events.PFJet.truth))
    l_mask = (truth == 0)

    if (iov == '2017') or (iov == '2018'):
        b_eff = 0.8  # overall btag eff for 2017/18
        light_eff = 0.03 # overall eff for light to be mistagged as b for 2017/18
    elif (iov == '2016') or (iov == '2016APV'):
        b_eff = 0.7 # overall btag eff for 2016
        light_eff = 0.04 # overall eff for light to be mistagged as b for 2016
    
    eff = ak.where(l_mask, light_eff, b_eff)
    
    # SF for VETO
    # denom
    denom_nottag = 1-(eff*np.ones(len(ak.flatten(mask_isTagged))))
    denom_nottag = ak.unflatten(denom_nottag, ak.num(events.PFJet.passMedID))

    denom_tag = (eff)*np.ones(len(ak.flatten(mask_isTagged)))
    denom_tag = ak.unflatten(denom_tag, ak.num(events.PFJet.passMedID))
    
    #denom = ak.where(mask_isTagged, 1, denom_nottag)
    denom = ak.where(mask_isTagged, denom_tag, denom_nottag)

    # num
    sf = get_btagSF(jsonPath, events.PFJet.pt, events.PFJet.eta, events.PFJet.truth, type)

    eff = ak.unflatten(eff, ak.num(events.PFJet.truth))
    
    num_tag = sf*eff
    #num_tag = ak.where(num_tag > 0, num_tag, eff)

    num_nottag = 1-sf*eff
    num_nottag = ak.where(num_nottag > 0, num_nottag, 1-eff)

    #num = ak.where(mask_isTagged, 1, num_nottag)
    num = ak.where(mask_isTagged, num_tag, num_nottag)

    #print('NUMERATOR 1-SF*eff or 1-eff based on sign', num)
    #print('DENOMINATOR 1-eff', denom)

    # SF for Tagging -> 1
    denom_final = ak.where(mask_isTagged, 1, denom_nottag)
    #print('denom_final', denom_final)
    denom_final_2 = ak.where(mask_isTagged, denom_tag, denom_nottag)

    num_final = ak.where(mask_isTagged, 1, num)
    num_final_2 = ak.where(mask_isTagged, num_tag, num_nottag)

    num_final = ak.prod(num_final, axis=-1)
    denom_final = ak.prod(denom_final, axis=-1)

    num_final_2 = ak.prod(num_final_2, axis=-1)
    denom_final_2 = ak.prod(denom_final_2, axis=-1)

    #print('final btag SF', num_final/denom_final)
    #print('final btag SF (accoutning for tag SF)', num_final_2/denom_final_2)

    return num_final/denom_final
    #return num_final_2/denom_final_2

    
##########
## JEC
###########

def applyJetCorrections(isMC, IOV, era=None, corr_type=None, do_factorized_jec_unc=False):
    jer_tag=None
    if (IOV=='2018'):
        jec_tag="Summer19UL18_V5_MC"
        jec_tag_data={
            "RunA": "Summer19UL18_RunA_V5_DATA",
            "RunB": "Summer19UL18_RunB_V5_DATA",
            "RunC": "Summer19UL18_RunC_V5_DATA",
            "RunD": "Summer19UL18_RunD_V5_DATA",
        }
        jer_tag = "Summer19UL18_JRV2_MC"
    elif (IOV=='2017'):
        jec_tag="Summer19UL17_V5_MC"
        jec_tag_data={
            "RunB": "Summer19UL17_RunB_V5_DATA",
            "RunC": "Summer19UL17_RunC_V5_DATA",
            "RunD": "Summer19UL17_RunD_V5_DATA",
            "RunE": "Summer19UL17_RunE_V5_DATA",
            "RunF": "Summer19UL17_RunF_V5_DATA",
        }
        jer_tag = "Summer19UL17_JRV2_MC"
    elif (IOV=='2016'):
        jec_tag="Summer19UL16_V7_MC"
        jec_tag_data={
            "RunF": "Summer19UL16_RunFGH_V7_DATA",
            "RunG": "Summer19UL16_RunFGH_V7_DATA",
            "RunH": "Summer19UL16_RunFGH_V7_DATA",
        }
        jer_tag = "Summer20UL16_JRV3_MC"
    elif (IOV=='2016APV'):
        jec_tag="Summer19UL16APV_V7_MC"
        ## HIPM/APV     : B_ver1, B_ver2, C, D, E, F
        ## non HIPM/APV : F, G, H
        jec_tag_data={
            "RunB-ver1": "Summer19UL16APV_RunBCD_V7_DATA",
            "RunB-ver2": "Summer19UL16APV_RunBCD_V7_DATA",
            "RunC": "Summer19UL16APV_RunBCD_V7_DATA",
            "RunD": "Summer19UL16APV_RunBCD_V7_DATA",
            "RunE": "Summer19UL16APV_RunEF_V7_DATA",
            "RunF": "Summer19UL16APV_RunEF_V7_DATA",
        }
        jer_tag = "Summer20UL16APV_JRV3_MC"
    else:
        raise ValueError(f"Error: Unknown year \"{IOV}\".")
    
    extract = extractor()
    corrPath = "../../../analysisTools/corrections/"
    if (isMC):
        #For MC
        extract.add_weight_sets([
            '* * ../../../analysisTools/corrections//JEC/{0}/{0}_L1FastJet_AK4PFchs.txt'.format(jec_tag),
            '* * ../../../analysisTools/corrections//JEC/{0}/{0}_L2Relative_AK4PFchs.txt'.format(jec_tag),
            '* * ../../../analysisTools/corrections//JEC/{0}/{0}_L3Absolute_AK4PFchs.txt'.format(jec_tag),
            '* * ../../../analysisTools/corrections//JEC/{0}/{0}_UncertaintySources_AK4PFchs.junc.txt'.format(jec_tag),
            '* * ../../../analysisTools/corrections//JEC/{0}/{0}_Uncertainty_AK4PFchs.junc.txt'.format(jec_tag),
        ])

        if jer_tag:
            extract.add_weight_sets([
            '* * ../../../analysisTools/corrections//JER/{0}/{0}_PtResolution_AK4PFchs.jr.txt'.format(jer_tag),
            '* * ../../../analysisTools/corrections//JER/{0}/{0}_SF_AK4PFchs.jersf.txt'.format(jer_tag)])
    else:       
        #For data, make sure we don't duplicate
        tags_done = []
        for run, tag in jec_tag_data.items():
            if not (tag in tags_done):
                extract.add_weight_sets([
                '* * ../../../analysisTools/corrections//JEC/{0}/{0}_L1FastJet_AK4PFchs.txt'.format(tag),
                '* * ../../../analysisTools/corrections//JEC/{0}/{0}_L2Relative_AK4PFchs.txt'.format(tag),
                '* * ../../../analysisTools/corrections//JEC/{0}/{0}_L3Absolute_AK4PFchs.txt'.format(tag),
                '* * ../../../analysisTools/corrections//JEC/{0}/{0}_L2L3Residual_AK4PFchs.txt'.format(tag),
                ])
                tags_done += [tag]
    
    extract.finalize()
    evaluator = extract.make_evaluator()
    
    if (isMC):
        jec_names = [
            '{0}_L1FastJet_AK4PFchs'.format(jec_tag),
            '{0}_L2Relative_AK4PFchs'.format(jec_tag),
            '{0}_L3Absolute_AK4PFchs'.format(jec_tag),
            '{0}_Uncertainty_AK4PFchs'.format(jec_tag)]
        if do_factorized_jec_unc:
            for name in dir(evaluator):
               #factorized sources
               if '{0}_UncertaintySources_AK4PFchs'.format(jec_tag) in name:
                    jec_names.append(name)
        if jer_tag: 
            jec_names.extend(['{0}_PtResolution_AK4PFchs'.format(jer_tag),
                              '{0}_SF_AK4PFchs'.format(jer_tag)])

    else:
        jec_names={}
        for run, tag in jec_tag_data.items():
            jec_names[run] = [
                '{0}_L1FastJet_AK4PFchs'.format(tag),
                '{0}_L3Absolute_AK4PFchs'.format(tag),
                '{0}_L2Relative_AK4PFchs'.format(tag),
                '{0}_L2L3Residual_AK4PFchs'.format(tag),]
    if isMC:
        jec_inputs = {name: evaluator[name] for name in jec_names}
    else:
        jec_inputs = {name: evaluator[name] for name in jec_names[era]}
    
    jec_stack = JECStack(jec_inputs)
    name_map = jec_stack.blank_name_map
    name_map['JetPt'] = 'pt'
    name_map['JetEta'] = 'eta'
    name_map['JetPhi'] = 'phi'
    name_map['JetMass'] = 'mass'
    name_map['Rho'] = 'rho'
    name_map['JetA'] = 'area'
    name_map['ptGenJet'] = 'pt_gen'
    name_map['ptRaw'] = 'pt_raw'
    name_map['massRaw'] = 'mass_raw'
    name_map['METpt'] = 'pt'
    name_map['METphi'] = 'phi'
    name_map['UnClusteredEnergyDeltaX'] = 'MetUnclustEnUpDeltaX'
    name_map['UnClusteredEnergyDeltaY'] = 'MetUnclustEnUpDeltaY'
    if corr_type=='met': return CorrectedMETFactory(name_map)
    return CorrectedJetsFactory(name_map, jec_stack)


def apply_JEC(isMC, IOV, events, type):
    # PFJet is already corrected, so uncorrect it first, and then re-correct it
    events.__setitem__("uncorrJet",events.PFJet)
    events["uncorrJet", "rho"]  = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, events.PFJet.pt)[0]
    events["uncorrJet", "pt_raw"] = (events.uncorrJet.rawFactor) * events.uncorrJet.pt
    events["uncorrJet", "mass_raw"] = (events.uncorrJet.rawFactor) * events.uncorrJet.mass
    events["uncorrJet", "pt"] = events.uncorrJet.pt_raw
    events["uncorrJet", "mass"] = events.uncorrJet.mass_raw
    events["uncorrJet", "pt_gen"] = events.uncorrJet.matchedGenJetPt
    events["uncorrJet", "p4"] = ak.with_name(events.uncorrJet[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")

    corrJetFactory = applyJetCorrections(isMC=True, IOV='2018', do_factorized_jec_unc=False)
    
    # overwrite PFJet with the uncorrected then recorrected version
    corrJet = corrJetFactory.build(events.uncorrJet, lazy_cache= events.caches[0])

    #print(corrJet.fields)
    if type == 'nominal':
        events.__setitem__("PFJet", corrJet)
    elif type == 'jes_up':
        events.__setitem__("PFJet", corrJet.JES_jes.up)
        correct_MET(events, type)
    elif type == 'jes_down':
        events.__setitem__("PFJet", corrJet.JES_jes.down)
        correct_MET(events, type)
    elif type == 'jer_up':
        events.__setitem__("PFJet", corrJet.JER.up)
        correct_MET(events, type)
    elif type == 'jer_down':
        events.__setitem__("PFJet", corrJet.JER.down)
        correct_MET(events, type)

def correct_MET(events, type):
    # PFJet is already corrected, so uncorrect it first, and then re-correct it
    print('PFMET pt before corr = ', events.PFMET.pt)
    if type == 'jes_up':
       events["PFMET","pt"] = events.PFMET.JESUpPt
       events["PFMET","phi"] = events.PFMET.JESUpPhi
    elif type == 'jes_down':
       events["PFMET","pt"] = events.PFMET.JESDownPt
       events["PFMET","phi"] = events.PFMET.JESDownPhi
    elif type == 'jer_up':
       events["PFMET","pt"] = events.PFMET.JERUpPt
       events["PFMET","phi"] = events.PFMET.JERUpPhi
    elif type == 'jer_down':
       events["PFMET","pt"] = events.PFMET.JERDownPt
       events["PFMET","phi"] = events.PFMET.JERDownPhi
    elif type == 'unclustered_up':
       events["PFMET","pt"] = events.PFMET.UnclusteredUpPt
       events["PFMET","phi"] = events.PFMET.UnclusteredUpPhi
    elif type == 'unclustered_down':
       events["PFMET","pt"] = events.PFMET.UnclusteredDownPt
       events["PFMET","phi"] = events.PFMET.UnclusteredDownPhi

    print('PFMET pt after corr = ', events.PFMET.pt)