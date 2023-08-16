import uproot
import sys

flist = sys.argv[1]

with open(flist,"r") as fin:
    files = fin.read().splitlines()

num_ev = 0
num_ev_tot = 0
for f in files:
    with uproot.open(f) as rf:
        num_ev += rf['Events'].num_entries
        num_ev_tot += 5000
print(float(num_ev)/float(num_ev_tot))