import numpy as np
import awkward as ak

def cut0(events,info):
    # using the medium WP, as in Andre's version of iDM
    name = "cut0"
    desc = "Dummy"
    plots = True
    return events, name, desc, plots