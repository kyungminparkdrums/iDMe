import os
import sys
import json
from tqdm import tqdm

input_json = sys.argv[1]

with open(input_json,"r") as fin:
    samples = json.load(fin)

output = {}

for sample in tqdm(list(samples.keys())):
    subsamples = samples[sample]
    for subsample in tqdm(list(subsamples.keys()),leave=False):
        dataset = subsamples[subsample]
        files = json.loads(os.popen(f"dasgoclient --query='file dataset={dataset}' -json").read())
        xrd = "root://cmsxrootd.fnal.gov/"
        files_out = {}
        for f in files:
            subF = f['file']
            for sf in subF:
                files_out[xrd+sf['name']] = sf['nevents']
        output[f"{sample}_{subsample}"] = files_out

fname = input_json.split("/")[-1].split(".")[0]
with open(f"{fname}_fileList.json","w") as fout:
    json.dump(output,fout,indent=4)