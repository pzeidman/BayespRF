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
3. Open the 'Run_first_level' script. Set the data_root_dir variable to the directory where you would like to store the example dataset. Then run the script. The GLM will be estimated and timeseries extracted for each of 10 runs.

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

The structure U will have one entry per time point (t) of the stimulus. E.g. if a bar crosses the display in 20 steps, and this happens 100 times, there will be 2000 entries in U. The U structure contains fields which describe the stimulus at time t:
- Dist and angle are vectors of the polar coordinates illuminated at this time step. 
- Ons and dur specify the onset and duration of this stimulus in seconds. Here we've set this to assume one frame of the stimulus per TR (the time taken to acquire an MRI volume), but this need not be the case.
- Modelling is conducted at a finer timescale than MRI acquisitions. Each stimulus time step will be divided into short 'microtime' bins. dt is the length of each bin in seconds (typically nmicrotime=16).
- pmax is the stimulus diameter in degrees of visual angle
- pmin is the minimum entertained pRF size in degrees of visual angle

*Note: The model has been tested with a stimulus resolution downsampled to a [41 x 41] grid. We recommend you do the same. See the example in scripts/prepare_inputs_polar_samsrf.m*

**2. Prepare your timeseries**

Using SPM to extract timeseries for pRF analysis is the recommended method. When viewing a single subject's SPM result, click Eigenvariate in the small grey window and follow the prompts. This will create a file named VOI_xx.mat, containing a summary timeseries, as well as timeseries for every voxel in the ROI. Create a cell array containing the filenames for each session's VOI:

```Matlab
xY = cell(1,num_sess);
for i = 1:num_sess
    filename = sprintf('VOI_Mask_%d.mat',sess(i));
    xY{i}    = fullfile(glm_dir,filename);
end
```

**3. Specify your pRF model**

Armed with your experiment timings U and timeseries xY, you can now specify your pRF model. This consists of creating an options structure with various settings, and calling the spm_prf_analyse function in 'specify' mode:

```Matlab
% Set pRF specification options
options = struct('TE', TE,...
                 'voxel_wise', true,...
                 'name', 'SamSrf_example',...
                 'model', 'spm_prf_fcn_gaussian_1sigma_DCP2',...
                 'B0',3);
             
% Specify pRF model (.mat file will be stored in the GLM directory)
PRF = spm_prf_analyse('specify',SPM,xY,U,options);
```
For a full list of options, see the help in spm_prf_analyse. Here we have just set a few key options: 
- TE is the echo time of our scanner, and is used in the neurosvascular model. 
- voxel_wise tells the model to work on a voxel-by-voxel basis, rather than using the summary (ROI) timeseries
- name will be used for the PRF file's filename
- model is the pRF model to use. There are a few to choose from (in the response_functions folder). Here we're using a circular pRF with polar coordinates
- B0 is the strength of the scanner's magnetic field in teslas. This influences several priors in the observation model.

This will create a file called PRF_name.mat in the same folder as the SPM.mat for this subject.

**4. Estimate the pRF model**

The evidence (free energy) and parameters of the model will now be estimated. For this example, we will only estimated one voxel:

```Matlab
% Estimation options
options  = struct('voxels',3439);

% Estimate
PRF_est = spm_prf_analyse('estimate',prf_file,options);
```

For a full set of estimation options, including the use of parallel computing, see the help in the spm_prf_analyse function.

**5. Review the results**

To review a pRF model's results, use the spm_prf_review function:

```Matlab
prf_file = fullfile(glm_dir,'PRF_SamSrf_example.mat');

spm_prf_review(prf_file, 3439);
```
Here we ask it to review only the voxel we have estimated (3439). If we don't give this praameter, and if the pRF file contains multiple voxel, then spm_prf_review will attempt to build a parameteric map of the voxels across the brain.
