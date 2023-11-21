# Condor jobs for AODSkimmer (`miniPlusNtuplizer`)
Condor jobs for running AODSkimmer on samples stored in FNAL LPC area.

## 1. Make filelists
Make the list of files containing the root file location in FNAL LPC area. In `../fileLists/`, 

```
python3 makeSignalFileLists.py [eos directory] 

# example
python3 makeSignalFileLists.py /store/group/lpcmetx/iDMe/Samples/signal/2018/AOD/
```

This will go over the subdirectories of the eos directory given as the argument. Then, it will make a `.txt` file of root file locations for each subdirectory, i.e. `Mchi-10p5_dMchi-1p0.txt` that contains the list of root files stored in `/store/group/lpcmetx/iDMe/Samples/signal/2018/AOD/Mchi-10p5_dMchi-1p0/`. 

You can also find the previously generated list of files in `AODSkimmer/fileLists/signal/`.

## 2. Set the output directory
In `miniPlusNtuplizer`, change in `exec_miniPlusNtuplizer.sh` (as necessary) the directory where output ntuples will be saved in [L23](https://github.com/kyungminparkdrums/iDMe/blob/main/AODSkimmer/condor/miniPlusNtuplizer/exec_miniPlusNtuplizer.sh#L23).

## 3. Submit condor job for a given filelist in `.txt`
In `miniPlusNtuplizer`, Run the following:
```
source submit_miniPlusNtuplizer_condor.sh [filelist in txt from Step 1] [year] [nThreads] [isData; 0/1] [isSignal; 0/1]

# example
source submit_miniPlusNtuplizer_condor.sh ../fileLists/signal/Mchi-10p5_dMchi-1p0_ctau-1.txt 4 0 1
```

## Some notes
- This will run the ntuplizer based on `ntuplizer_CMSSW_10_6_26.tar.gz` in `/store/group/lpcmetx/iDMe/compiled_CMSSW_envs/` at the LPC eos area. If you have made any changes to your ntuplizer, make sure that you extract your `CMSSW` environment in `.tar.gz` and that your `CMSSW` environment is used in condor jobs. This can be done by updating [L11](https://github.com/kyungminparkdrums/iDMe/blob/main/AODSkimmer/condor/miniPlusNtuplizer/exec_miniPlusNtuplizer.sh#L11C1-L11C6) in `exec_miniPlusNtuplizer.sh`.
- In `miniPlusNtuplizer/miniPlusNtuplizer_config.jdl`, you will need to change the directory of `executable`, `output`, `error`, and `log` to match yours. 
