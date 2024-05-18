from CRABAPI.RawCommand import crabCommand
import os
import sys
from io import StringIO

directory = sys.argv[1]

original_stdout = sys.stdout

sys.stdout = StringIO()
s = crabCommand("status",dir=directory)
tmp_stdout = sys.stdout
sys.stdout = original_stdout
outlines = tmp_stdout.getvalue().splitlines()
link = None
for line in outlines:
    if "Dashboard monitoring URL:" in line:
        for i in range(len(line)):
            if line[i:i+5] == "https":
                link = line[i:]
                break
    if link is not None:
        break
print(directory)
print(f"\t STATUS : {s['status']}")
for k in s['jobsPerStatus'].keys():
    print(f"\t {k} : {s['jobsPerStatus'][k]}")
print(f"\t LINK : {link}")
print("--------------------------")