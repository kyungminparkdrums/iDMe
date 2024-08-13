import os
import glob

process_list = glob.glob('crab*')

for process in process_list:
    if '.log' in process:
        continue
    elif 'DYJetsToLL_M-4to50' in process:
        continue
    else:
        cmd = 'crab resubmit -d {} --maxmemory=2500 --maxjobruntime=2750'.format(process)

        print(cmd)
        os.system(cmd)

