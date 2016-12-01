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

%% Prepare inputs

% Build a structure containing which stimulus pixels were illuminated at
% each time step.
load(fullfile(data_dir,'aps_Bars.mat'));
U = prepare_inputs_polar_samsrf(ApFrm,TR);

% The timeseries from each session are stored in a VOI_xx.mat file. Build a
% cell array of the VOI files for each session.
xY = cell(1,num_sess);
for i = 1:num_sess
    filename = sprintf('VOI_Mask_%d.mat',sess(i));
    xY{i}    = fullfile(glm_dir,filename);
end
%% Specify pRF model (all 6669 voxels)

% Load SPM for timing information / image dimensions
SPM = load(fullfile(glm_dir,'SPM.mat'));
SPM = SPM.SPM;

% Update SPM path as we don't know where this example will be saved
SPM.swd = glm_dir;

% Set pRF specification options
options = struct('TE', TE,...
                 'voxel_wise', true,...
                 'name', 'SamSrf_example',...
                 'model', 'spm_prf_fcn_gaussian_1sigma_DCP2',...
                 'B0',3);
             
% Specify pRF model (.mat file will be stored in the GLM directory)
PRF = spm_prf_analyse('specify',SPM,xY,U,options);

%% Estimate one voxel (global peak) as an example: idx 3439, [-4,-70,-3] 

% Model to estimate
prf_file = fullfile(glm_dir,'PRF_SamSrf_example.mat');

% Estimation options
options  = struct('voxels',3439);

% Estimate
PRF_est = spm_prf_analyse('estimate',prf_file,options);

%% Review (single voxel)
prf_file = fullfile(glm_dir,'PRF_SamSrf_example.mat');

spm_prf_review(prf_file, 3439);

%% Merge
P = {};
for i = 1:186
    P{i} = sprintf('models_wb/PRF_SamSrf_example_%d.mat',i);
end
spm_prf_analyse('merge',P,pwd);