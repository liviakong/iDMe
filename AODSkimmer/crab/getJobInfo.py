from CRABAPI.RawCommand import crabCommand
import os
import sys
from io import StringIO

directory = sys.argv[1]
dirs = [d for d in os.listdir(directory) if os.path.isdir(f"{directory}/{d}")]

original_stdout = sys.stdout

for d in dirs:
    sys.stdout = StringIO()
    s = crabCommand("status",dir=f"{directory}/{d}")
    sys.stdout = original_stdout
    print(d)
    print(f"\t STATUS : {s['status']}")
    for k in s['jobsPerStatus'].keys():
        print(f"\t {k} : {s['jobsPerStatus'][k]}")
    print("--------------------------")