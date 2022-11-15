#!/bin/bash


cmsRun scripts/run_ElectronNtuplizer_cfg.py year=2018 flist="file:mini_Mchi-48p0_dMchi-16p0_ctau-1.root" outfile="regIDNoCut_Mchi-48p0_dMchi-16p0_ctau-1.root"
cmsRun scripts/run_ElectronNtuplizer_cfg.py year=2018 flist="file:mini_Mchi-48p0_dMchi-16p0_ctau-10.root" outfile="regIDNoCut_Mchi-48p0_dMchi-16p0_ctau-10.root"
cmsRun scripts/run_ElectronNtuplizer_cfg.py year=2018 flist="file:mini_Mchi-48p0_dMchi-16p0_ctau-100.root" outfile="regIDNoCut_Mchi-48p0_dMchi-16p0_ctau-100.root"
cmsRun scripts/run_ElectronNtuplizer_cfg.py year=2018 flist="file:mini_Mchi-48p0_dMchi-16p0_ctau-1000.root" outfile="regIDNoCut_Mchi-48p0_dMchi-16p0_ctau-1000.root"

hadd regIDNoCut_Mchi-48p0_dMchi-16p0_combined.root regIDNoCut_Mchi-48p0_dMchi-16p0_ctau-1.root regIDNoCut_Mchi-48p0_dMchi-16p0_ctau-10.root regIDNoCut_Mchi-48p0_dMchi-16p0_ctau-100.root regIDNoCut_Mchi-48p0_dMchi-16p0_ctau-1000.root
