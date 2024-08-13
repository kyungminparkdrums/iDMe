import numpy as np
import awkward as ak
import sys

'''
def cut1(events,info):
    name = "cut1"
    desc = "Preselections (no HEM veto) & n(good_vtx) & n(Jet) Cuts"
    plots = True

    return events, name, desc, plots
'''

def cut2(events,info):
    name = "cut2"
    desc = "HEM Jet Veto (eta upper bound bug)"
    plots = True

    if len(events) != 0:
        cut = events.hasHEMjetBug == 0
        print(f'HEM Jet Veto Pass (with eta bug): {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []

    return events[cut], name, desc, plots

def cut3(events,info):
    name = "cut3"
    desc = "HEM Jet Veto (additionally veto missing eta region)"
    plots = True

    if len(events) != 0:
        cut = events.hasHEMjet == 0
        print(f'HEM Jet Veto Pass (additional eta region): {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []

    return events[cut], name, desc, plots

def cut4(events,info):
    name = "cut4"
    desc = "HEM PF electron Veto"
    plots = True

    if len(events) != 0:
        cut = (events.hasHEMelecPF == 0)
        print(f'HEM PF electron Veto Pass: {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []

    return events[cut], name, desc, plots

def cut5(events,info):
    name = "cut5"
    desc = "HEM Lpt electron Veto"
    plots = True

    if len(events) != 0:
        cut = (events.hasHEMelecLpt == 0)
        print(f'HEM Lpt electron Veto Pass: {np.count_nonzero(cut)}/{len(cut)}')
    else:
        cut = []

    return events[cut], name, desc, plots
