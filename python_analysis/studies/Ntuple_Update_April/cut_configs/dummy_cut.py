import numpy as np
import awkward as ak
import sys
sys.path.append("../analysisTools/")
import analysisSubroutines as routines

model_dir = "./models/"

def cut0(events,info):
    # pre-computing BDT score so I can make plots with it
    name = 'cut4'
    desc = 'dummy'
    plots = True
    
    if len(events) != 0:
        input = routines.makeBDTinputs(events)

        model = f'{model_dir}/BDT_inclusive_10Vars_minDxyCut.json'
        score_BDT = routines.getBDTscore(input, model)

        events['BDTScore'] = score_BDT
    else:
        events['BDTScore'] = ak.zeros_like(events.eventWgt)

    return events, name, desc, plots