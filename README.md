# BayespRF Toolbox
The BayespRF Toolbox is a framework for specifying and comparing models of population receptive fields (pRFs). Models are probabilistic and are estimated using the variational Bayes algorithm. The toolbox is built in Matlab and depends upon [SPM](http://www.fil.ion.ucl.ac.uk/spm/).

Advantages:

- Neural and neuro-vascular / BOLD parameters are estimated on a per-voxel (or region of interest) basis
- Models for any stimulus dimensionality and with any number of parameters can be estimated
- The uncertainty (variance) of parameters is estimated as well as the covariance among parameters
- Models expressing competing hypotheses can be compared based on their evidence (using the free energy approximation).

Disadvantages:

- Model estimation is computationally expensive (~100 seconds per voxel per CPU core)
- Cortical surface projection is not yet implemented

*The toolbox is very new and still in development. So please enjoy using it, but treat results with caution*

## Installing
1. Copy the 'toolbox' folder to a location on your computer of your choice.
2. Add the 'toolbox' folder to your Matlab path with subfolders.

## Running the example dataset (retinotopic mapping)
A good way to get started is to try fitting a pRF model using an example dataset. The toolbox includes a folder called **example_3T** which contains the necessary scripts to automatically download and analyse the example dataset. The dataset is a single subject scanned with 3T fMRI, kindly shared online by the [SamSrf](https://figshare.com/articles/SamSrf_toolbox_for_pRF_mapping/1344765) toolbox.

**Part 1: Downloading the data, extracting timeseries**

1. Copy the example_3T folder, included with BayespRF, to a location of your choice on your computer.
2. Navigate to the **example_3T/scripts** folder in Matlab.
3. Open the **Run_first_level** script. Set the data_root_dir variable to the directory where you would like to store the example dataset. Then run the script. The GLM will be estimated and timeseries extracted for each of 10 runs.

**Part 2: Running the pRF analysis**

1. Edit the **example_3T/scripts/Run_pRF_analysis** script. Change the **data_root_dir** variable to match the location containing the example dataset.
2. Run the Run_pRF_analysis script. A pRF model will be specified covering 6669 voxels. A single pRF will then be estimated (index 3439) and a window will be displayed with the results.

A step-by-step walkthrough of the demo can be found below.

## What's in the toolbox

The toolbox provides a set of functions for specifying and analysing pRF models:

| Function | Description |
| -------- | ----------- |
| spm_prf_analyse | Creates and estimates pRF models. |
| spm_prf_bpa_within_subject | Performs Bayeian Parameter Averaging (BPA) at the within-subject level |
| spm_prf_editor | A graphical editor for adjusting priors in pRF models |
| spm_prf_find_voxel | Gets the index of a voxel within a pRF based on its mm coordinates |
| spm_prf_get_ppd | Gets the prior and posterior predictive density for a pRF model |
| spm_prf_review | Reviews the results of pRF estimation |

A pRF model is defined by a response function - for example, a Gaussian response or Difference of Gaussians (DoG) response. pRF models are provided in the **toolbox/response_functions** folder. The name of the function to use should be provided when specifying the pRF in **spm_prf_analyse**.

### Response functions

| | Function | Description | Input coordinates |
| --- | -------  | ----------- | ----------------- |
| ![Gaussian](https://cloud.githubusercontent.com/assets/2145293/20843434/41822ca2-b8b3-11e6-8ccc-e473dffb648f.png) | spm_prf_fcn_gaussian_1sigma_DCP2 | Isotropic 2D Gaussian | Polar |
| ![Ellipitical Gaussian](https://cloud.githubusercontent.com/assets/2145293/20843432/4181f37c-b8b3-11e6-8438-ee4346e23bea.png) | spm_prf_fcn_gaussian_1sigma_ellipse_DCP2 | Elliptical 2D Gaussian | Polar |
| ![Rotated Gaussian](https://cloud.githubusercontent.com/assets/2145293/20843433/4181f872-b8b3-11e6-94e1-7175f5d44a27.png) | spm_prf_fcn_gaussian_1sigma_angled_DCP2 | Elliptical 2D Gaussian with rotation | Polar |
| ![DoG](https://cloud.githubusercontent.com/assets/2145293/20843436/4182ab0a-b8b3-11e6-818a-94d2dfd45fe9.png) | spm_prf_fcn_gaussian_DoG_basic_DCP2 | Isotropic 2D DoG | Polar |
| ![Ellipitical DoG](https://cloud.githubusercontent.com/assets/2145293/20843437/418757c2-b8b3-11e6-9975-9dba3c6284be.png) | spm_prf_fcn_gaussian_DoG_ellipse_DCP2 | Elliptical 2D DoG | Polar |
| ![Rotated DoG](https://cloud.githubusercontent.com/assets/2145293/20843438/4192a348-b8b3-11e6-80ac-06745b46a49b.png) | spm_prf_fcn_gaussian_DoG_DCP2 | Elliptical 2D DoG with rotation | Polar |

The neurovascular signal model (**spm_prf_fx.m**) and the BOLD signal model (**spm_prf_gx.m**) do not need to be modified on a study-by-study basis.

## Step by Step

Here's a step by step guide to running a pRF analysis. We'll be working through the steps in the example script supplied with the toolbox (scripts/Run_pRF_analysis.m).

**1. Prepare your inputs**

You need to give BayespRF the timing of your experimental stimuli - in other words, which pixels of the stimuli were illuminated at which time points. For the polar coordinates model used in this example, this will be a Matlab structure with the following fields:

```Matlab
    U(t).dist  = dist;         % Distance
    U(t).angle = angle;        % Angle
    U(t).ons = TR * (t-1);     % Onset (secs)
    U(t).dur = TR;             % Duration (secs)
    U(t).dt  = TR/nmicrotime;  % 1 bin=1/16 second, i.e. 16 bins per second
    U(t).pmax = stim_diameter; % Stimulus diameter
    U(t).pmin = 0.5;           % Minimum PRF size
```

The structure U will have one entry per time point (t) of the stimulus. Time points are indexed from 1. For example, if a bar crosses the display in 20 steps, and this happens 100 times, there will be 2000 entries in U. 

The U structure contains fields which describe the stimulus at time t:
- **Dist** and **angle** are vectors of the polar coordinates illuminated at this time step. 
- **Ons** and **dur** specify the onset and duration of this stimulus in seconds. Here we've set this to assume one frame of the stimulus per TR (the time taken to acquire an MRI volume), but this need not be the case.
- **dt** is the length of each microtime bin in seconds (typically nmicrotime=16). Modelling is conducted at a finer timescale than MRI acquisitions. Each stimulus time step will be divided into short 'microtime' bins of length dt seconds.
- **pmax** is the stimulus diameter in degrees of visual angle
- **pmin** is the minimum entertained pRF size in degrees of visual angle

**Note \#1:** You will find an example script for producing this U structure in **scripts/prepare_inputs_polar_samsrf.m** . This reads in the stimuli as a 3D matrix of stimuli, where the third dimension indexes time, and produces the U structure.  

**Note \#2:** The pRF models have been tested with the resolution of the stimuli downsampled to a [41 x 41] grid. We recommend you do the same with your stimuli. This is done for you in the example script.

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

Armed with your experiment timings U and timeseries xY, you can now specify your pRF model. This consists of creating an options structure with various settings, and calling the **spm_prf_analyse** function in 'specify' mode:

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
- **TE** is the echo time of the scanning sequence, and is used in the neurosvascular model. 
- **voxel_wise** tells the model to work on a voxel-by-voxel basis, rather than using the summary (ROI) timeseries
- **name** will be used for the PRF's filename
- **model** is the pRF model to use. There are a few to choose from (in the response_functions folder). Here we're using a circular pRF with polar coordinates
- **B0** is the strength of the scanner's magnetic field in teslas. This influences several priors in the observation model.

This will create a file called PRF_name.mat in the same folder as the SPM.mat for this subject.

**4. Estimate the pRF model**

The evidence (free energy) and parameters of the model will now be estimated. For this example, we will only estimated one voxel:

```Matlab
% Estimation options
options  = struct('voxels',3439);

% Estimate
PRF_est = spm_prf_analyse('estimate',prf_file,options);
```

For a full set of estimation options, including the use of parallel computing, see the help in the **spm_prf_analyse** function.

**5. Review the results**

To review a pRF model's results, use the spm_prf_review function:

```Matlab
prf_file = fullfile(glm_dir,'PRF_SamSrf_example.mat');

spm_prf_review(prf_file, 3439);
```
Here we ask it to review only the voxel we have estimated (3439). If we don't give this parameter, and if the pRF file contains multiple voxels, then spm_prf_review will attempt to build a parameteric map of the voxels across the brain.

