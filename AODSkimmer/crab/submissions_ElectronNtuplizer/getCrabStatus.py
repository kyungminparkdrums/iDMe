from CRABAPI.RawCommand import crabCommand
import os
import glob
import datetime

directories = glob.glob('./crab_*/')

list_completed = []
list_failed = []
list_running = []

now = datetime.datetime.now()

f = open("./JobStatus_{}_{}_{}_{}_{}_{}.txt".format(now.year, now.month, now.day, now.hour, now.minute, now.second),"w")

for idx, directory in enumerate(directories):
    s = crabCommand("status",dir=directory)
    
    f.write("\n\nAnalyzing "+directory)
    f.write("\nstatus = " + s["status"])

    #process = directory.split('crab_iDMe_')[-1].split('_')[0]

    if s['status'] == "COMPLETED":
        list_completed.append(directory)
        f.write("\nJobs COMPLETED!")
        os.system("mv {} ./SUCCEEDED_JOBS".format(directory))
    elif s['status'] == "FAILED":
        list_failed.append(directory)
        failedJobs = [j[1] for j in s['jobList'] if j[0] == 'failed']
        f.write("\nJobs FAILED! Failed job IDs = ")
        for j in failedJobs:
            f.write(j + ", ")
    else:
        list_running.append(directory)
        f.write("\nJobs still running!")

f.close()


flist = open("./ProcessLists_{}_{}_{}_{}_{}_{}.txt".format(now.year, now.month, now.day, now.hour, now.minute, now.second),"w")

flist.write("\n\nProcess with crab job completed: \n")
for p in list_completed:
    flist.write(p + "\n")

flist.write("\n\nProcess with crab job failed: \n")
for p in list_failed:
    flist.write(p + "\n")

flist.write("\n\nProcess with crab job still running: \n")
for p in list_running:
    flist.write(p + "\n")

flist.close()

