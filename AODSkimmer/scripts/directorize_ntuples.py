from XRootD import client

prefix = "/store/group/lpcmetx/iDMe/Samples/Ntuples/"
xrdClient = client.FileSystem("root://cmseos.fnal.gov")
status, flist = xrdClient.dirlist(prefix)
rootFiles = [item.name for item in flist if '.root' in item.name]

for f in rootFiles:
    name = f.split(".")[0]
    mchi = ""
    dmchi = ""
    ct = ""
    
    for s in name.split("_"):
        if "Mchi" in s and "dMchi" not in s:
            mchi = s.split("-")[1]
        elif "dMchi" in s:
            dmchi = s.split("-")[1]
        elif "ctau" in s:
            ct = s.split("-")[1]
       
    xrdClient.mkdir(prefix+"signal/2018/Mchi-{0}_dMchi-{1}/ctau-{2}/".format(mchi,dmchi,ct),flags=client.flags.MkDirFlags.MAKEPATH)
    target = prefix+"signal/2018/Mchi-{0}_dMchi-{1}/ctau-{2}/{3}".format(mchi,dmchi,ct,f)
    xrdClient.mv(prefix+f,target)
