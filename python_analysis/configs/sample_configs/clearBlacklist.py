import json
import sys

infile = sys.argv[1]
with open(infile,"r") as f:
    samples = json.load(f)
for s in samples:
    s['blacklist'] = []
with open(infile,"w") as f:
    json.dump(samples,f,indent=4)