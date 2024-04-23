import os
from coffea.processor import accumulate
from coffea.util import save, load
import hist
import sys
from XRootD import client
import shutil
import numpy as np
from tqdm import tqdm

def copy_histos(loc,floc="tmp_merge"):
    if os.path.isdir(f"./{floc}/"):
        shutil.rmtree(f"./{floc}/")
    os.mkdir(f"{floc}")
    xrdClient = client.FileSystem("root://cmseos.fnal.gov")
    status, flist = xrdClient.dirlist(loc)
    files = ["root://cmseos.fnal.gov/"+loc+"/"+f.name for f in flist if ".coffea" in f.name]
    for f in files:
        os.system(f"xrdcp {f} ./{floc}/")    

def merge_histos(loc,sample_type,floc="tmp_merge",n_seg=3):
    files = os.listdir(floc)
    
    outDirName = f"{sample_type}/{loc.split('/')[-1]}/" if loc[-1] != "/" else f"{sample_type}/{loc.split('/')[-2]}/"
    print(outDirName)
    if os.path.isdir(outDirName):
        shutil.rmtree(outDirName)
    os.mkdir(outDirName)
    
    """print("getting list of histograms")
    hkeys = load(f"tmp_merge/{files[0]}")[0]
    out_keys = list(hkeys.keys())
    del hkeys

    for i,k in enumerate(out_keys):
        print(f"Merging {k}, {i+1}/{len(out_keys)}")
        toMerge = []
        for f in files:
            hin = load(f"tmp_merge/{f}")[0]
            toMerge.append(hin[k].copy())
            del hin
        outObj = accumulate(toMerge)
        save(outObj,f"{outDirName}/{k}.coffea")
        del outObj, toMerge"""
        

    """histograms = []
    metadata = []
    files = os.listdir("tmp_merge/")
    print("loading initial file")
    histos, meta = load(f'tmp_merge/{files[0]}')
    ifile = 1
    for f in files[1:]:
        print(f"loading file {ifile}")
        h_in, meta_in = load(f'tmp_merge/{f}')
        print(f'merging file {ifile}')
        histos = accumulate([histos,h_in])
        print(f'merging metadata {ifile}')
        meta = accumulate([meta,meta_in])
        print("deleting loaded files")
        del h_in, meta_in
        ifile += 1
    save((histos,meta),f"{outDirName}/histograms.coffea")"""

    test = load(f"{floc}/{files[0]}")[0]
    hnames = list(test.keys())
    del test
    segments = np.array_split(files,n_seg)
    for j,seg in enumerate(segments):
        print("segment",j+1,"of",len(segments))
        histograms = []
        print('loading files')
        for f in tqdm(seg):
            histograms.append(load(f"{floc}/{f}")[0])
        print('accumulating histos')
        #for i,k in enumerate(tqdm(hnames)):
        #    histos = accumulate([h[k] for h in histograms])
        #    save(histos,f"{outDirName}/{k}_seg{j}.coffea")
        #    del histos
        hout = accumulate(histograms)
        save(hout,f"{outDirName}/seg_{j}.coffea")
        del hout
        del histograms
    
    """print('merging segments')
    for i,k in enumerate(tqdm(hnames)):
        hout = []
        for j in range(len(segments)):
            hout.append(load(f"{outDirName}/{k}_seg{j}.coffea"))
            os.remove(f"{outDirName}/{k}_seg{j}.coffea")
        hout = accumulate(hout)
        save(hout,f"{outDirName}/{k}.coffea")
        del hout
    print('merged histos')"""
    

if __name__ == "__main__":
    loc = sys.argv[1]
    sample_type = sys.argv[2]
    copy = bool(int(sys.argv[3]))
    merge = bool(int(sys.argv[4]))
    n_seg = int(sys.argv[5])
    if len(sys.argv) > 6:
        floc = sys.argv[6]
    print("copy",copy)
    print("merge",merge)
    if copy:
        print("copying")
        copy_histos(loc,floc=floc)
    if merge:
        merge_histos(loc,sample_type,floc=floc,n_seg=n_seg)