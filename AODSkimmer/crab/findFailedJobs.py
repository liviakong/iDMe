from CRABAPI.RawCommand import crabCommand
import os
import sys

directory = sys.argv[1]

print "Analyzing "+directory
s = crabCommand("status",dir=directory)
if s['status'] != "COMPLETED" and s['status'] != "FAILED":
    print "status = "+s["status"]
    print "Jobs still running!"
    exit
failedJobs = [j[1] for j in s['jobList'] if j[0] == 'failed']
f = open(directory+"/failedJobs.txt","w")
for j in failedJobs:
    f.write(j+"\n")
f.close()