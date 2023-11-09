# Crab job for AODSkimmer

## 1. Set up environment
```
source setupAPI.sh
cmsenv
``` 

## 2. Submit crab jobs
```
python submit_multicrab_miniPlusNtuplizer.py -c submit -y [year] -f [file list in json] -t [isData 0/1] -s [isSignal 0/1]

# example
python submit_multicrab_miniPlusNtuplizer.py -c submit -y 2018 -f ../../fileLists/background/bkg_QCD_2018.json -t 0 -s 0
```

