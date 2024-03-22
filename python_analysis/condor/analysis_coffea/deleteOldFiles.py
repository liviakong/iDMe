from XRootD import client
import sys
from datetime import datetime as dt

pref = sys.argv[1]
num_days = float(sys.argv[2])

xrdClient = client.FileSystem("root://cmseos.fnal.gov")
files = [f for f in xrdClient.dirlist(pref,flags=client.flags.DirListFlags.STAT)[1] if '.coffea' in f.name]
now = dt.now()
old = [f for f in files if (now-dt.fromisoformat(f.statinfo.modtimestr)).total_seconds() >= num_days*86400]
dts = [(now-dt.fromisoformat(f.statinfo.modtimestr)).total_seconds() for f in files]
names = [f.name for f in files]
for f in files:
    print(f.name)
    print(f.statinfo.modtimestr)
for f in old:
    print(pref+'/'+f.name)
    print(f.statinfo.modtimestr)
