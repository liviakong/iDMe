## Cut configs `cut_configs`
Cuts applied in coffea output are defined in a cut config in `.py`. [histo_configs/SR_v2_skimmed](https://github.com/kyungminparkdrums/iDMe/blob/main/python_analysis/configs/cut_configs/SR_v2_skimmed.py) is an example config file for running analyzer on skimmed (preselected) ntuples. Inside the cut config, each cut is defined in a function. Inside the function, you can set `plots` variable to `True` or `False` based on whether you want to save histograms after that cut or not. You can change the cut value in `cut` variable inside each function. Note that when you turn on saving the plots for all steps, sometimes coffea might throw out memory error (if you have "too many" histograms saved).

## Histo configs `histo_configs`

When running coffea analyzer, a histogram config in `.py` will be needed. [histo_configs/SR_studies.py](https://github.com/kyungminparkdrums/iDMe/blob/main/python_analysis/configs/histo_configs/SR_studies.py) is an example config file. 

Histograms saved in the coffea output are filled through [fillHistos](https://github.com/kyungminparkdrums/iDMe/blob/main/python_analysis/configs/histo_configs/SR_studies.py#L107) function in the config file, which is called in the [analysis tools](https://github.com/kyungminparkdrums/iDMe/blob/iDMe/python_analysis/analysisTools/analysisTools.py#L133) when running coffea. 

Binning of the histograms are taken care of in [histobins.py](https://github.com/kyungminparkdrums/iDMe/blob/main/python_analysis/configs/histo_configs/histobins.py), which is imported in histo config files.

### Add histograms in coffea output
If you want to add any histograms in the coffea output file, edit the following:

1. In `histo_configs/histobins.py`, add the variable and its range & binning. 
```
# example
dR = Regular(100,0,6,name='dr',label="$\Delta R$")
vtx_matchType = IntCategory([0,1,2],name="mtype",label="Vertex Gen Match Type")
``` 

2. In a histo config json file you want to use, i.e. `histo_configs/SR_studies.py`, add the histogram in `make_histograms` function inside `histograms` dictionary.
```
# example
"sel_vtx_dR" : Hist(samp,cut,dR,storage=hist.storage.Weight()),                                     # 1D
"sel_vtx_minEledRj_vs_matchType" : Hist(samp,cut,dR,vtx_matchType,storage=hist.storage.Weight()),   # 2D
```

Inside the argument for constructing `Hist` instance, the third (for 1D; for 2D, it'll be third and 4th) argument after `cut` should match the variable you added in `histobins.py`. 

3. In the histo config file, add in `fillHistos` function what the histograms should be filled with. 
```
# example
histos["sel_vtx_dR"].fill(samp=samp,cut=cut,dr=vtx.dR,weight=wgt), 
histos["sel_vtx_minEledRj_vs_matchType"].fill(samp=samp,cut=cut,dr=np.minimum(e1.mindRj,e2.mindRj),mtype=vtx.match,weight=wgt),
```

Inside the argument for `fill`, th third (for 1D; for 2D, it'll be third and 4th) argument after `cut` should match the `name` of the variable you added in `histobins.py`. The value you assign to the `name` of the variable is named after the branch name in the ntuples. For example, in the ntuples we have [vtx_dR](https://github.com/kyungminparkdrums/iDMe/blob/iDMe/CustomTools/src/NtupleContainerV2.cc#L207C8-L207C8) branch. Coffea will read it as `vtx.dR`, and we assign this to `dr` in this case.  


## Sample configs `sample_configs`

### Create json config for signal samples

1. Make a json file with signal samples in eos area 
In `sample_configs`, 
```
python makeSignalConfigs.py [mode: sig or bkg] [year] [alpha] [prefix for ntuple location] [output json file name]

# example
python makeSignalConfigs.py sig 2018 aEM /store/group/lpcmetx/iDMe//Samples/Ntuples/signal_v2/ test_signal_v2
```

The example command will create `test_signal_v2_2018_aEM.json` that contains the location of each signal samples stored in `/store/group/lpcmetx/iDMe//Samples/Ntuples/signal_v2/` subdirectories. But `xsec`, `sum_wgt`, `blacklist` will be empty.

2. Fill out the cross section information in the json file.
In `sample_configs`,
```
python getSignalXsec.py [json file from Step 1]

# example
python getSignalXsec.py test_signal_v2_2018_aEM.json
```

3. Fill out the sum of gen weights information and blacklist 'bad' files.
In `sample_configs`,
```
python sumGenWgts.py [json file from Step 2]

# example
python sumGenWgts.py test_signal_v2_2018_aEM.json
```

