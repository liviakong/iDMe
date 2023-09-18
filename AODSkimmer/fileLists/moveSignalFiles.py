from XRootD import client
from XRootD.client.flags import MkDirFlags as mkdflags
import re

cmseos  = "root://cmseos.fnal.gov"
xrd = "root://cmsxrootd.fnal.gov/"
xrdClient = client.FileSystem(cmseos)

base = "/store/group/lpcmetx/iDMe//Samples/signal/2018/AOD/"
files = [item.name for item in xrdClient.dirlist(base)[1] if '.root' in item.name]

for f in files:
    m_dm = re.search("Mchi-(\d+p\d)_dMchi-(\d+p\d)",f)[0]
    ctau = int(float(re.findall("ctau-(\d+\.\d)",f)[0]))
    ctau = f'ctau-{ctau}'
    destination = f"{base}/{m_dm}/{ctau}/"
    status, items = xrdClient.dirlist(destination)
    if status.error: # directory doesn't exist yet
        print(f"Making directory {destination}")
        xrdClient.mkdir(destination,flags=mkdflags.MAKEPATH)
    xrdClient.mv(f"{base}/{f}",f"{destination}/{f}")
