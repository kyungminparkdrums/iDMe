import os
import glob

flist = glob.glob('../../../fileLists/signal/2018/MINIAOD/*.txt')

for f in flist:
    cmd = f'source ./submit_ElectronNtuplizer_condor.sh {f} 2018 4 0 1 v12_miniAOD 1'
    print(cmd)
    os.system(cmd)
