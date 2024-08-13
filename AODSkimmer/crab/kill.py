import os
import glob

process_list = glob.glob('crab*')

for process in process_list:
    if '.log' in process:
        continue
    else:
        cmd = 'crab kill -d {}'.format(process)

        print(cmd)
        os.system(cmd)

