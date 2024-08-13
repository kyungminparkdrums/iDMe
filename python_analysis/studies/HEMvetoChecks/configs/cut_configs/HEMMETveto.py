import numpy as np
import awkward as ak
import sys


def cut4(events,info):
    # HEM MET veto
    name = "cut4"
    desc = "HEM MET veto"
    plots = True
    hasHEMMET = (events.PFMET.phi < -0.87) & (events.PFMET.phi > -1.57)
    # (met.pt>470)|(met.phi>-0.62)|(met.phi<-1.62)
    cut = hasHEMMET == 0
    return events[cut], name, desc, plots
