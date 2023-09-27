# Coffea analysis tools for iDMe

This directory contains all the python-based analysis tools for iDMe. The analysis pipeline (i.e. applying cuts, filling histograms, etc.) is built using the [coffea](https://coffeateam.github.io/coffea/) framework for columnar analysis. Here's a summary of everything in this area:
- `analysisTools` : core functionality for the `coffea`-based analysis, including code for the `Analyzer` and `Processor` objects.
- `configs` : houses three different kinds of config files we use with our coffea processors: histograms, cuts, and samples. The dedicated README inside this directory explains things in more detail
- `condor` : tools for submitting `coffea` analysis jobs to condor. See README inside for details
- `studies` : houses directories and notebooks for various studies. Feel free to add your own stuff here!
- `testNotebooks` :  