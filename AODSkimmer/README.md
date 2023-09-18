# AODSkimmer

This folder contains all the code to produce ntuples from AOD. The core analysis script is `plugins/ElectronSkimmer.cc`, which contains all the code necessary for our analysis (primarily recording event information (e.g. MET) and vertexing pairs of electrons). One quirk of our analysis is that we run on top of AOD, rather than miniAOD. This is due to problems with the default miniAOD event content for low-$p_T$ electrons, which does not include the correct energy regression and uses a sub-optimal version of the BDT ID.

## Running the analysis code locally
The script for running the analysis code is `scripts/miniPlusElectronNtuplizer_cfg.py`. Here's an example of how to run it (make sure you run `cmsenv` first!):

```bash
cmsRun scripts/miniPlusElectronNtuplizer_cfg.py data=0/1 signal=0/1 \ 
year=YEAR numThreads=N \
nEvents=1000 flist=PATH_TO_FILE_LIST outfile=my_output.root
```
The `data=0/1` argument specifies whether or not you're running on data (should be `0` for all Monte Carlo samples, whether signal or background). Similarly, `signal=0/1` specifies if you're running over a signal sample or not. For example: `data=0 signal=1` would be for running over signal MC, `data=0 signal=0` would be for background MC, and `data=1 signal=0` would be for data. The `year` argument specifies what data/MC year the sample you're running over corresponds to, e.g. `year=2018` for running on 2018 data samples or MC samples with the 2018 configuration. The `numThreads` argument specifies how many threads to use for multithreading (defaults to 8). The `nEvents` argument specifies how many events to process (defaults to `-1`, in which case all events from the input files will be processed). The `flist` argument should be a path to a `.txt` file with a list of ROOT file locations for whatever sample you're analyzing. Lastly, the `outfile=my_output.root` specifies a name for the output root file (defaults to `test_output.root`).

## Example command
Here's an example command to analyze 1000 events from the $m_\chi = 10$ GeV, $\Delta m_\chi = 1$ GeV, $c\tau = 1$ mm signal sample
```bash
cmsRun scripts/miniPlusElectronNtuplizer_cfg.py year=2018 data=0 signal=1 \
nEvents=1000 flist="fileLists/signal/Mchi-10p5_dMchi-1p0_ctau-1.txt" \
outfile=test_signal_output.root
```

## File lists
File lists for signal and background samples are found in sub-folders of the `fileLists` directory. Currently, it is only possible to run over signal files on the FNAL LPC (either locally or using condor). Running over backgrounds requires submitting jobs to `crab` (see the README in the `crab` directory).

## Condor job submission
See the README in the `condor` directory for instructions on how to submit analysis jobs to condor.