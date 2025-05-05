import subprocess
import glob
import os

indirs = glob.glob(f'bkg_*/')
#indirs = glob.glob(f'data_*/')

for indir in indirs:
    #print(indir)
    subdirs = glob.glob(f'{indir}/*/')
    for sub in subdirs:
        try:
            with open(f'{sub}/Logs/stdout.out', 'r') as f:
                lines = [line.rstrip() for line in f]
            final = lines[-1]

            if 'filtered' not in final:
            #if 'filtered' in final:
                #print(final)
            #else:
                #print(final)
                print(f'no filtered in {sub}')
        except:
            print(f'no stdout.out in {sub}')

