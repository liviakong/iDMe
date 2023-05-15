from XRootD import client
from XRootD.client.flags import MkDirFlags as mkdflags
import re
import sys

if len(sys.argv) < 2:
    print("Error: need to specify base directory")
    exit()

base = sys.argv[1]
cmseos  = "root://cmseos.fnal.gov"
xrd = "root://cmsxrootd.fnal.gov/"
xrdClient = client.FileSystem(cmseos)
mass_pts = [item.name for item in xrdClient.dirlist(base)[1] if '.root' not in item.name]

for mass_pt in mass_pts:
    lifetimes = [item.name for item in xrdClient.dirlist(f"{base}/{mass_pt}/")[1] if '.root' not in item.name]
    for ct in lifetimes:
        out_file = f"{mass_pt}_{ct}.txt"
        root_files = [f"{xrd}{base}/{mass_pt}/{ct}/{item.name}" for item in xrdClient.dirlist(f"{base}/{mass_pt}/{ct}/")[1] if '.root' in item.name]
        with open(out_file,"w") as fout:
            for f in root_files:
                fout.write(f+"\n")


"""
#!/bin/bash

basedir=$1

cmseos="root://cmseos.fnal.gov/"
xrd="root://cmsxrootd.fnal.gov/"

for d in `eos $cmseos ls $basedir`
do
    filename="${d}.txt"
    touch $filename
    for f in `eos $cmseos ls $basedir/$d/*.root`
    do
        echo $xrd$basedir/$d/$f >> $filename
    done
done
"""