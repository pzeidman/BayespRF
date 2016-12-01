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

## Step by Step

Here's a step by step guide to running a pRF analysis. We'll be working through the steps in the example script supplied with the toolbox (scripts/Run_pRF_analysis.m).

**1. Prepare your inputs**

You need to give BayespRF the timing of your experimental stimuli - in other words, which pixels of the stimuli were illuminated at which time points. For the polar coordinates model used in this example, this will be a Matlab structure as follows:

```Matlab
    U(t).dist  = dist;         % Distance
    U(t).angle = angle;        % Angle
    U(t).ons = TR * (t-1);     % Onset (secs)
    U(t).dur = TR;             % Duration (secs)
    U(t).dt  = TR/nmicrotime;  % 1 bin=1/16 second, i.e. 16 bins per second
    U(t).pmax = stim_diameter; % Stimulus diameter
    U(t).pmin = 0.5;           % Minimum PRF size
```

The structure U will have one entry per time point (t) of the stimulus. E.g. if a bar crosses the display in 20 steps, and this happens 100 times, there will be 2000 entries in U. The structure contains fields which describe the stimulus at time t:
- Dist and angle are vectors of the polar coordinates illuminated at this time step. 
- Ons and dur specify the onset and duration of this stimulus in seconds. Here we've set this to assume one frame of the stimulus per TR (the time taken to acquire an MRI volume), but this need not be the case.
- Modelling is conducted at a finer timescale than MRI acquisitions. Each stimulus time step will be divided into short 'microtime' bins. dt is the length of each bin in seconds (typically nmicrotime=16).
- pmax is the stimulus diameter in degrees of visual angle
- pmin is the minimum entertained pRF size in degrees of visual angle
