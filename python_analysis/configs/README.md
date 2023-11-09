# Sample configs

## Create json config for signal samples

### 1. Make a json file with signal samples in eos area 
In `sample_configs`, 
```
python makeSignalConfigs.py [mode: sig or bkg] [year] [alpha] [prefix for ntuple location] [output json file name]

# example
python makeSignalConfigs.py sig 2018 aEM /store/group/lpcmetx/iDMe//Samples/Ntuples/signal_v2/ test_signal_v2
```

The example command will create `test_signal_v2_2018_aEM.json` that contains the location of each signal samples stored in `/store/group/lpcmetx/iDMe//Samples/Ntuples/signal_v2/` subdirectories. But `xsec`, `sum_wgt`, `blacklist` will be empty.

### 2. Fill out the cross section information in the json file.
In `sample_configs`,
```
python getSignalXsec.py [json file from Step 1]

# example
python getSignalXsec.py test_signal_v2_2018_aEM.json
```

### 3. Fill out the sum of gen weights information and blacklist 'bad' files.
In `sample_configs`,
```
python sumGenWgts.py [json file from Step 2]

# example
python sumGenWgts.py test_signal_v2_2018_aEM.json
```

