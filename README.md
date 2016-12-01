# BayespRF
The BayespRF Toolbox is a framework for specifying and comparing models of population receptive fields (pRFs). It is built in Matlab and depends upon [SPM](http://www.fil.ion.ucl.ac.uk/spm/).

## Installing
1. Copy the 'toolbox' folder to a location on your computer of your choice.
2. Add the 'toolbox' folder to your Matlab path with subfolders.

## Running the example dataset
A good way to get started is to try fitting a pRF model using an example dataset. The toolbox includes a folder called 'example_3T' which contains the necessary scripts to automatically download and analyse the example dataset supplied with the [SamSrf](https://figshare.com/articles/SamSrf_toolbox_for_pRF_mapping/1344765) toolbox.

**Part 1: Downloading the data, extracting timeseries**

1. Copy the 'example_3T' folder, included with BayespRF, to a location of your choice on your computer.
2. Navigate to the 'example_3T/scripts' folder in Matlab.
3. Edit the 'Run_first_level' script - set the data_root_dir variable to the directory where you would like to store the example dataset.
4. Run the 'Run_first_level' script. The GLM will be estimated and timeseries extracted for each of 10 runs.

**Part 2: Running the pRF analysis**

1. Edit the Run_pRF_analysis script, which you will find in the 'example_3T/scripts' folder. Change the data_root_dir variable to match the location containing the example dataset.
2. Run the Run_pRF_analysis script. A pRF model will be specified covering 6669 voxels. A single pRF will then be estimated (index 3439) and a window will be displayed with the results.
