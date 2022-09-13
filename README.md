# cceps_single_trials

This repository contains all of the code needed to reproduce the figures from Cornblath et al. 2022 ("Quantifying trial-by-trial variability during cortico-cortical evoked potential mapping of epileptic networks"). Data can be found on ieeg.org.
 
## General overview

There are two main components of this pipeline:

	- `cceps_files.m` and `cceps_files.R`: these files allow you to specify your own custom paths
	- `trial_variability/main.m` : this script runs all scripts necessary to process CCEP data and calculate N1 and N2 for each stimulation trial
	- `trial_variability/pub_figs/main.R` : this script runs all the main analyses of CCEP data for the manuscript.

Please contact Eli Cornblath (Eli ~`DOT`~ Cornblath ~`AT`~ pennmedicine.upenn.edu) with any questions regarding this code.