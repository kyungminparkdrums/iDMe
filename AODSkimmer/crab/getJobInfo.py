from CRABAPI.RawCommand import crabCommand
import os
import sys
from io import StringIO

directory = sys.argv[1]

original_stdout = sys.stdout

sys.stdout = StringIO()
s = crabCommand("status",dir=directory)
tmp_stdout = sys.stdout
sys.stdout = original_stdout
outlines = tmp_stdout.getvalue().splitlines()
link = None
for line in outlines:
    if "Dashboard monitoring URL:" in line:
        for i in range(len(line)):
            if line[i:i+5] == "https":
                link = line[i:]
                break
    if link is not None:
        break
print(directory)
print(f"\t STATUS : {s['status']}")
for k in s['jobsPerStatus'].keys():
    print(f"\t {k} : {s['jobsPerStatus'][k]}")
print(f"\t LINK : {link}")
failed_all = []
failed_main = []
failed_tail = []
for j in s['jobList']:
    status,name = j[0], j[1]
    if status == 'failed':
        if name[:2] == '0-':
            continue
        failed_all.append(name)
        if name[:2] == '1-':
            failed_tail.append(name)
        if "-" not in name:
            failed_main.append(name)
true_failed = []
for f in failed_main:
    if f"1-{f}" in failed_tail:
        true_failed.append(f)
for f in failed_tail:
    if f[2:] not in failed_main:
        true_failed.append(f)
print("Failed jobs:",failed_all)
print("Failed main:",failed_main)
print("Failed tail:",failed_tail)
print("True failed:",",".join(true_failed))
print(f"crab resubmit {directory} --jobids={','.join(true_failed)}")
print("--------------------------")
