
## Condor jobs for skimmer that applies preselections to the ntuples

We use AOD samples for the analysis, so we run the [ntuplizer+miniAOD] from [AODSkimmer](https://github.com/kyungminparkdrums/iDMe/tree/main/AODSkimmer). We call these *unskimmed ntuples*. One could run the coffea analyzer on the unskimmed ntuples (using the [sample json files](https://github.com/kyungminparkdrums/iDMe/tree/main/python_analysis/configs#sample-configs). 

One could also make *skimmed ntuples* where some [basic preselections](https://github.com/kyungminparkdrums/iDMe/blob/main/python_analysis/condor/condor_skim_rdf.py#L182-L188) are applied to the *unskimmed ntuples*.  

### 1. Change (as necessary) output directory of the skimmer.
In [L32](https://github.com/kyungminparkdrums/iDMe/blob/main/python_analysis/condor/submit_condor_skim_rdf.py#L32), change the output directory of the skimmed ntuples.

### 2. Submit condor jobs for the skimmer.
```
python submit_condor_skim_rdf.py [sample json file] [n_files_per] [n_cores]

# example
python submit_condor_skim_rdf.py ../configs/sample_configs/signal_v2_2018_aEM.json 10 4
```

The example command will go over the signal files in `../configs/sample_configs/signal_v2_2018_aEM.json`, and make an output root file processing over each 10 input root files. For example, if there are 20 unskimmed input files, 2 skimmed output files will be generated. 

## Condor jobs for running coffea analyzer on ntuples (either *unskimmed ntuples* or *skimmed ntuples*)
