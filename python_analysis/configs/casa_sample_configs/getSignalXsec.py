import pandas as pd
import sys
import json

inputJson = sys.argv[1]

df = pd.read_csv('/uscms/home/sbrightt/nobackup/iDMe/signal_xsec/MG5_aMC_v2_9_6/bin/signal_xsec_table.csv')
with open(inputJson) as f:
    samples = json.load(f)
for samp in samples:
    mchi = samp["Mchi"]
    dmchi = samp["dMchi"]
    ct = samp["ctau"]
    aD = str(samp['alphaD'])
    sel_row = df[(df["Mchi"] == mchi) & (df["dMchi"] == dmchi) & (df["ct"] == ct) & (df["alphaD"] == aD)]
    if sel_row.empty:
        print(f"No xsec found for {samp['name']}")
        samp["xsec"] = 0.0
    else:
        samp["xsec"] = float(sel_row["xsec(pb)"].iloc[0]) * 1000 # convert to fb

with open(inputJson,'w') as f:
    json.dump(samples,f,indent=4)