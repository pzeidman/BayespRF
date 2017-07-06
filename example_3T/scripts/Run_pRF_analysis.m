% Example script to run a simple pRF analysis on a voxel-by-voxel basis 
% and view the results.
%
% Neuronal model:      Single Gaussian function
% Receptive field:     Circular (isotropic)
% Input specification: Polar coordinates
%
% Data is from SamSrf

% Settings

% Directory of the downloaded example dataset
data_root_dir = 'D:\samdata\';
data_dir      = fullfile(data_root_dir,'example','pRF');

% Directory of GLM
glm_dir  = fullfile(pwd,'../GLM');

TR = 1;         % Repetition time
TE = 0.055;     % Echo time

% Which sessions to include
sess = 1:10;
num_sess = length(sess);

% The hemisphere to analyse. We'll just do left for now.
hemi = 'lh';

%% Prepare inputs

% Build a structure containing which stimulus pixels were illuminated at
% each time step.
load(fullfile(data_dir,'aps_Bars.mat'));
U = prepare_inputs_polar_samsrf(ApFrm,TR);

% The timeseries from each session are stored in a VOI_xx.mat file. Build a
% cell array of the VOI files for each session.
xY = cell(1,num_sess);
for i = 1:num_sess
    filename = sprintf('VOI_%s_prf_mask_%d.mat',hemi,sess(i));
    xY{i}    = fullfile(glm_dir,filename);
end
%% Specify pRF model (all 2422 voxels)

% Load SPM for timing information / image dimensions
SPM = load(fullfile(glm_dir,'SPM.mat'));
SPM = SPM.SPM;

% Update SPM path as we don't know where this example will be saved
SPM.swd = glm_dir;

% Set pRF specification options
options = struct('TE', TE,...
                 'voxel_wise', true,...
                 'name', [hemi '_SamSrf_example'],...
                 'model', 'spm_prf_fcn_gaussian_polar',...
                 'B0',3);
             
% Specify pRF model (.mat file will be stored in the GLM directory)
PRF = spm_prf_analyse('specify',SPM,xY,U,options);

%% Estimate one voxel as an example (voxel 1)
voxel = 1;

% Model to estimate
prf_file = fullfile(glm_dir,['PRF_' hemi '_SamSrf_example.mat']);

% Estimation options
options  = struct('voxels',voxel);

% Estimate
PRF_est = spm_prf_analyse('estimate',prf_file,options);

% Review
spm_prf_review(prf_file, voxel);
%% Estimate all voxels (slow)

% Model to estimate
prf_file = fullfile(glm_dir,['PRF_' hemi '_SamSrf_example.mat']);

% Estimation options
options  = struct('use_parfor',true);

% Estimate
PRF_est = spm_prf_analyse('estimate',prf_file,options);

% Review
spm_prf_review(prf_file);
%% Convert the left V1 label to a nifti mask so we can do some ROI analyses
label_file = fullfile(data_dir, 'lh_V1.label');
spm_prf_import_label( label_file, glm_dir );

%% Plot the summed pRF response in right V1

% Load estimated pRF file
prf_file = fullfile(glm_dir,['PRF_' hemi '_SamSrf_example.mat']);
load(prf_file);

% Load VOI (imported using spm_prf_import_label)
roi = fullfile(glm_dir, [hemi '_V1.nii']);

figure('Color','w');
spm_prf_summarise(PRF,roi);
title('Region of interest','FontSize',16);
%% Compute a negative entropy map (certainty of the pRF location)

% Load estimated pRF file
prf_file = fullfile(glm_dir,['PRF_' hemi '_SamSrf_example.mat']);
load(prf_file);

% Compute and plot
spm_prf_plot_entropy(PRF,{'dist','angle'},'dist_angle',true);