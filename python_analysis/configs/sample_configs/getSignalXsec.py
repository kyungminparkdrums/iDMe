import pandas as pd
import sys
import json

inputJson = sys.argv[1]
kind = sys.argv[2]
if kind != "sig" and kind != "bkg":
    print("Bad kind input ",kind)
    print("Use 'sig' or 'bkg'")

with open(inputJson) as f:
    samples = json.load(f)

if kind == 'sig':
    df = pd.read_csv('/uscms_data/d3/sbrightt/iDMe/signal_xsec/condor/signal_xsec_table.csv')
    for samp in samples:
        mchi = samp["Mchi"]
        dmchi = samp["dMchi"]
        ct = samp["ctau"]
        aD = str(samp['alphaD'])
        sel_row = df[(df["Mchi"] == mchi) & (df["dMchi"] == dmchi) & (df["ct"] == ct) & (df["alphaD"] == aD) & (df['mA/m1'] == 3)]
        if sel_row.empty:
            print(f"No xsec found for {samp['name']}")
            samp["xsec"] = 0.0
        else:
            samp["xsec"] = float(sel_row["xsec(pb)"].iloc[0]) * 1000 # convert to fb

if kind == 'bkg':
    with open("bkg_xsecs.json","r") as fin:
        xsecs = json.load(fin)
    for samp in samples:
        name = samp['name']
        if name not in list(xsecs.keys()):
            print(f"No xsec for sample {name}")
            samp["xsec"] = 0.0
        else:
            samp["xsec"] = xsecs[name]

with open(inputJson,'w') as f:
    json.dump(samples,f,indent=4)