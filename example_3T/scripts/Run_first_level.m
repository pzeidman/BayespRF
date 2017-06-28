%% Downloads example data (900MB), runs GLM analysis and extracts timeseries

% Settings

% Directory into which to download example dataset
data_root_dir = 'D:\samdata\';
data_dir      = fullfile(data_root_dir,'example','pRF');
surf_dir      = fullfile(data_root_dir,'example','surf');

% Directory for creating GLM
glm_dir  = fullfile(pwd,'../GLM');

% Number of sessions
nsess = 10;

% Repetition time
TR    = 1;
%% Download and unzip example data
if ~exist(data_root_dir,'file')
    mkdir(data_root_dir);
end

% Download
fn = fullfile(data_root_dir,'Example.zip');
fprintf('%-40s:', 'Downloading SamSrf dataset...');
urlwrite('https://zenodo.org/record/163582/files/Example.zip',fn);

% Unzip
unzip(fn, data_root_dir);

fprintf(' %30s\n', '...done');
%% Prepare onsets
load(fullfile(data_dir,'aps_Bars.mat'));
U = prepare_inputs_polar_samsrf(ApFrm,TR);

bins_x = [-8.5 0 8.5];
bins_y = [8.5 0 -8.5];

% Build time x screen bins matrix (3x3 screen bins)
onsets_matrix = zeros(length(U), length(bins_x) .^ 2);
for t = 1:length(U)
    % Loop over pixels activated at this time point
    for activated_pixel = 1:length(U(t).dist)
        % Get location
        dist  = U(t).dist(activated_pixel);
        angle = U(t).angle(activated_pixel);
        
        % Polar->cartesian
        x = dist * cos(angle);
        y = dist * sin(angle);

        % Identify closest bin
        [~,binx] = min(abs(bins_x-x));
        [~,biny] = min(abs(bins_y-y));

        % Binned coordintes -> index
        bin_idx = sub2ind([length(bins_x) length(bins_x)],biny,binx);

        onsets_matrix(t,bin_idx) = onsets_matrix(t,bin_idx) + 1;
    end
end

% Remove empty bins
onsets_matrix = onsets_matrix(:,any(onsets_matrix));
num_regressors = size(onsets_matrix,2);

% SPM inputs
names = cell(1,num_regressors); 
onsets = cell(1,num_regressors); 
durations = cell(1,num_regressors); 

for t = 1:num_regressors
    names{t} = ['Bin' num2str(t)];
    onsets{t} = (find( onsets_matrix(:,t) ) - 1) * TR;
    durations{t} = 0;
end

save('onsets.mat', 'names', 'onsets', 'durations');

%% Specify first level design

start_dir = pwd;

% Make output directory
if ~exist(glm_dir,'file')
    mkdir(glm_dir);
end

load('first_level_batch.mat');

% Session-specific options
for i = 1:nsess
    movement = spm_select('FPList',data_dir,sprintf('rp_bfBars%d.txt',i));
    epis     = spm_select('ExtFPList',data_dir,sprintf('ubfBars%d.nii',i), 1:999);    
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = cellstr(movement);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans     = cellstr(epis);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi     = cellstr('onsets.mat');
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;    
end

% Model spec options
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(glm_dir);
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;

% Run
spm_jobman('run',matlabbatch);

cd(start_dir);
%% Import cortical surface
%
% Creates images: GLM/lh_surface.nii and GLM/rh_surface.nii
% and .mat files: GLM/lh_Srf.mat and GLM/rh_Srf.mat

% Structural image
struct = fullfile(data_dir,'T1.nii');

% Left hemisphere
spm_prf_import_surface(glm_dir, struct, surf_dir, 'lh');

% Right hemisphere
spm_prf_import_surface(glm_dir, struct, surf_dir, 'rh');
%% Build a mask of voxels which survive p < 0.001
clear matlabbatch;
matlabbatch{1}.spm.stats.results.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.export{1}.binary.basename = 'mask_uncorrected';
spm_jobman('run',matlabbatch);
%% Remove voxels from the mask anterior to y = 0
cd(glm_dir);

% Read
V = spm_vol('spmF_0001_mask_uncorrected.nii');
[Y,XYZmm] = spm_read_vols(V);

% Threshold
i = XYZmm(2,:) > 0;

% Write
Y(i) = 0;
spm_write_vol(V,Y);

cd(start_dir);
%% Extract timeseries from surface voxels which survive p < 0.001

% Hemisphere (we'll just do left for now)
hemi = 'lh';

% Identify masks
spm_F_mask   = fullfile(glm_dir,'spmF_0001_mask_uncorrected.nii');
surface_mask = fullfile(glm_dir,[hemi '_surface.nii']);

% Prepare batch
load('extract_timeseries_batch.mat');
matlabbatch{1}.spm.util.voi.name   = [hemi '_prf_mask'];
matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(spm_F_mask);
matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr(surface_mask);
matlabbatch{1}.spm.util.voi.expression        = 'i1 & i2';

% Run batch
spm_jobman('run',matlabbatch);