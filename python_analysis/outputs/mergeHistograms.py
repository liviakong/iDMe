import os
from coffea.processor import accumulate
from coffea.util import save, load
import hist
import sys
from XRootD import client
import shutil

loc = sys.argv[1]
sample_type = sys.argv[2]

os.mkdir("tmp_merge")
xrdClient = client.FileSystem("root://cmseos.fnal.gov")
status, flist = xrdClient.dirlist(loc)
files = ["root://cmseos.fnal.gov/"+loc+"/"+f.name for f in flist if ".coffea" in f.name]
for f in files:
   os.system(f"xrdcp {f} ./tmp_merge/")

histograms = []
metadata = []
for f in os.listdir('tmp_merge/'):
    histos, meta = load(f'tmp_merge/{f}')
    histograms.append(histos)
    metadata.append(meta)
histograms = accumulate(histograms)
metadata = accumulate(metadata)

outFileName = f"{sample_type}/{loc.split('/')[-1]}.coffea"
save((histograms,metadata),outFileName)
shutil.rmtree("./tmp_merge/")



