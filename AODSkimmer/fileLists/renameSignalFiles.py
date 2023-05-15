from XRootD import client
from XRootD.client.flags import MkDirFlags as mkdflags
import re

cmseos  = "root://cmseos.fnal.gov"
xrd = "root://cmsxrootd.fnal.gov/"
xrdClient = client.FileSystem(cmseos)

base = "/store/group/lpcmetx/iDMe//Samples/signal/2018/AOD/"
files = [item.name for item in xrdClient.dirlist(base)[1] if '.root' in item.name]

for f in files:
    ct = re.findall("1.973269804e-(\d+).root",f)[0]
    print(f)
    if int(ct) == 13:
        new = f.replace("1.973269804e-13.root","1.0_year-2018.root")
        print(f"found {ct}, rename to {new}")
    elif int(ct) == 14:
        new = f.replace("1.973269804e-14.root","10.0_year-2018.root")
        print(f"found {ct}, rename to {new}")
    elif int(ct) == 15:
        new = f.replace("1.973269804e-15.root","100.0_year-2018.root")
        print(f"found {ct}, rename to {new}")
    elif int(ct) == 16:
        new = f.replace("1.973269804e-16.root","1000.0_year-2018.root")
        print(f"found {ct}, rename to {new}")
    #ctau = int(float(re.findall("ctau-(\d+\.\d)",f)[0]))
    #ctau = f'ctau-{ctau}'
    #destination = f"{base}/"
    #status, items = xrdClient.dirlist(destination)
    #if status.error: # directory doesn't exist yet
    #    print(f"Making directory {destination}")
    #    xrdClient.mkdir(destination,flags=mkdflags.MAKEPATH)
    xrdClient.mv(f"{base}/{f}",f"{base}/{new}")
