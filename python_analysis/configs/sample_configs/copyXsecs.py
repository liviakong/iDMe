import json
import sys

fsrc = sys.argv[1]
fref = sys.argv[2]

with open(fsrc,"r") as fin:
    dsrc = json.load(fin)
with open(fref,"r") as fin:
    dref = json.load(fin)

xsec_ref = {k['name']:k['xsec'] for k in dref}
for entry in dsrc:
    entry['xsec'] = xsec_ref[entry['name']]

with open(fsrc,"w") as fout:
    json.dump(dsrc,fout,indent=4)