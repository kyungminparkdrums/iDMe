import json
import os
import subprocess
with open('signal_2018.json') as f:
    filedict = json.load(f)

new_dict = {}
for samp in filedict.keys():
    new_dict[samp] = {}
    for ctau in [1,10,100]:
        point = f'{samp.split("sig_")[-1]}_ctau-{str(ctau)}'
        loc = filedict[samp][point]
        campaign = loc.split("/")[-2]
        wildcard = loc.replace(campaign,"*UL17*")
        cmd = f'dasgoclient --query="dataset={wildcard}"'
        output = subprocess.check_output(cmd, shell=True)
        newloc = output.decode("utf-8").split("\n")[0]
        print(newloc)
        new_dict[samp][point] = newloc 
    print(samp)

out_file = open('signal_2017.json', "w")
json.dump(new_dict, out_file, indent=4)
out_file.close()


#dasgoclient --query="dataset=/iDM_DarkPhotonToEE_Mchi-5p25_dMchi-0p5_ctau-100_mA-15p0_HT80_TuneCP5_13TeV_madgraph-pythia8/*UL17*/MINIAODSIM"
