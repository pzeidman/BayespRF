function varargout = spm_prf_analyse(mode,varargin)
% Runs a PRF analysis
%
% Examples and default settings
%
% -------------------------------------------------------------------------
% Specify an analysis:
%
% PRF = spm_prf_analysis('specify', SPM, xY, U, options);
%
% Where:
%
% SPM - filename of the subject's SPM.mat file (or SPM structure)
%
% xY  - Cell array of VOI_XX.mat files containing the timeseries for pRF
%       estimation, or a cell array of structures. Each element of the cell 
%       array is one session (run) with the format:
%
%       xY{s}.Y     - summary timeseries [t x 1]
%       xY{s}.y     - timeseries for each of v voxels [t x v]
%       xY{s}.XYZmm - locations of each timeseries y in the brain [3 x v]
%
%       For a voxel-wise analysis (options.voxel_wise == true), the data in
%       xY{s} will be used. For an ROI analysis (options.voxel_wise ==
%       false) the single timeseries from xY.Y will be used, which by
%       default is the principal eigenvariate of the ROI.
%
%       For multiple sessions, if options.avg_sess == true, the timeseries
%       of each session will be averaged. If optoins.avg_sess == false, the
%       timeseries will be concatenated.
%
% U   - Experimental stimulus timing. This is a struct array of size
%       [1 x u] for u stimulus timepoints. The fields are as follows:
%
%       U(t).dist     - polar distance coordinate of each pixel p
%                       illuminated at time t. Vector of size [p x 1].
%       U(t).angle    - polar angle coordinate of each pixel p
%                       illuminated at time t. Vector of size [p x 1].
%
%       U(t).pmax     - diameter of the stimulated area in degress of
%                       visual angle [1 x 1]
%       U(t).pmin     - minimum allowed pRF width (standard deviation) in
%                       degrees of visual angle [1 x 1]
%
%       U(t).onset    - stimulus onset time in seconds [1 x 1]
%       U(t).duration - stimulus duration in seconds [1 x 1]
%
%       U(t).dt       - each stimulus time step will be divided into short
%                       'microtime' bins. This is the length of each bin in 
%                       seconds, typically 1/8. [1 x 1]
%
%       Note that dist, angle, pmax and pmin may be redefined by the pRF
%       response function.
%
% options - Structure of pRF model specification options. E.g.
%           options = ...
%              struct('TE',TE, ...        % echo time
%                  'voxel_wise',false, ...% per voxel (true) or ROI (false)
%                  'model', model,...     % pRF function (spm_prf_fcn...)
%                  'hE',6, ...            % expected log precision of noise
%                  'P',[], ...            % starting parameters
%                  'B0',3, ...            % fMRI field strength (teslas)
%                  'avg_sess',true, ...   % average sessions' timeseries
%                  'avg_method','mean'    % accepts 'mean' or 'eigen'
%                  'delays', TR / 2 ...   % microtime offset
%              );
%
% -------------------------------------------------------------------------
% Estimate a pRF model:
%
% PRF = spm_prf_analyse('estimate', PRF, options);
%
% PRF      - Specified PRF model, containing one or more timeseries
% options  - Structure of estimation options. E.g.
%
%            options = ...
%               struct('init', 'GLM_P', ...     % Initialization. See below
%                      'use_parfor', false, ... % Parallelization
%                      'nograph', false, ...    % If true, disables plots
%                      'voxels', 1:v ....       % Voxels indices (optional)
%               );
%
%            The init parameter controls whether to initialize the
%            parameters using a simple grid search. This speeds up the
%            subsequent Bayesian estimation. Options:
%
%            'GLM_P' - Initializes the Bayesian estimation based on the
%            grid search (default).
%
%            'GLM' - Sets the priors of the Bayesian estimation based on a
%            grid search.
%
%            'NONE' - Initialize parameters at their default prior values.
%
%            The use_parfor function requires the Matlab Computing Toolbox.
%            Voxels' timeseries will be estimated in separate threads.
%
% -------------------------------------------------------------------------
% Split a multi-voxel pRF model into separate files containing one or more
% voxels each:
%
% spm_prf_analyse('split', PRF, nmodels, dir_out);
%
% PRF     - the specified pRF model
% nmodels - the target number of files to split the model into
% dir_out - output directory for the split models
%
% This function is particularly useful in a cluster setting, for producing
% separate packets of work for individual nodes to perform.
%
% -------------------------------------------------------------------------
% Merge pRF model files containing one or more voxels each into a single
% multi-voxel pRF file:
%
% spm_prf_analyse('merge', P, dir_out);
%
% P - cell or character array of pRF files to merge
% dir_out - output directory for the merged model
%
% This function is particularly useful in a cluster setting, for merging 
% pRFs estimated on different nodes in a cluster.
%                       
% -------------------------------------------------------------------------
% Extract a single voxel pRF from a multi-voxel pRF file:
%
% PRF = spm_prf_analyse('extract',PRF,voxel_idx,name);
%
% PRF       - filename or pRF structure
% voxel_idx - index of the voxel to extract
% name      - saves the output as PRF_name.mat
%
% ---------------------------------------------------------------------
% Copyright (C) 2016 Peter Zeidman
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.   
% ---------------------------------------------------------------------

switch upper(mode)
    case 'SPECIFY'
        SPM = varargin{1};
        xY  = varargin{2};
        U   = varargin{3};
        if length(varargin) >= 4
            options = varargin{4};
        else
            options = struct();
        end
        PRF = specify(SPM,xY,U,options);
        varargout{1} = PRF;
    case 'ESTIMATE'
        PRF     = varargin{1};
        try
            options = varargin{2};
        catch
            options = struct();
        end
        PRF = estimate(PRF,options);
        varargout{1} = PRF;
    case 'SPLIT'
        PRF        = varargin{1};
        num_models = varargin{2};
        if length(varargin) > 2
            out_dir = varargin{3};
        else
            out_dir = '';
        end
        
        if ischar(PRF)
            if isempty(out_dir), out_dir = fileparts(PRF); end
            PRF = load(PRF);
            PRF = PRF.PRF;                        
        elseif isempty(out_dir)
            out_dir = pwd;
        end
        split_model(PRF,num_models,out_dir);
    case 'MERGE'
        P = varargin{1};
        if ~iscell(P)
            P = cellstr(P);
        end        
        if length(varargin) > 1
            out_dir = varargin{2};
        else
            try out_dir = fileparts(P{1}); catch, out_dir = pwd; end
        end
        merge_models(P,out_dir);
    case 'EXTRACT'
        PRF      = varargin{1};
        idx      = varargin{2};
        out_name = varargin{3};
        try out_dir = fileparts(PRF); catch, out_dir = pwd; end
        
        if ischar(PRF)
            PRF = load(PRF);
            PRF = PRF.PRF;
        end
        
        PRF      = extract_model(PRF,idx);
        PRF.name = ['PRF_' out_name];
        
        save_prf(PRF, out_dir);
        varargout{1} = PRF;
    case 'GET_PARAMETERS'
        PRF = varargin{1};
        
        varargout{1} = get_corrected_parameters(PRF);
    otherwise
        error('Unknown action');
end
% -------------------------------------------------------------------------
function PRF = specify(SPM,rois,U,options)
% Specify a PRF model
%
% SPM     - SPM.mat, for getting the TR and working directory
% rois    - filename or structure with eigenvariate Y and xY (see 
%           spm_regions). If a cell array is given, the timeseries from 
%           each are averaged.
% U       - input timing
% options - model options

Y = [];
XYZmm = [];

% Set defaults
try options.voxel_wise; catch, options.voxel_wise = false; end
try options.hE;         catch, options.hE = 6; end
try options.P;          catch, options.P = []; end
try options.B0;         catch, options.B0 = 3; end

% Load SPM
if ischar(SPM)
    SPM = load(SPM);
    SPM = SPM.SPM;
end

if ~iscell(rois)
    rois = {rois};
end

n_sess = length(rois);

% Choose whether to average timeseries
try 
    options.avg_sess;   
catch
    options.avg_sess = (n_sess > 1);
end

% Choose method of averaging multiple timeseries
try 
    options.avg_method
catch
    options.avg_method = 'mean';
end

for sess = 1:n_sess
    
    % Unpack
    if ischar(rois{sess})
        % From filename
        voi         = load(rois{sess});
        xY_sess     = voi.xY;
        Y_sess      = voi.Y;
    elseif isstruct(rois{sess}) && isfield(rois{sess},'Y')
        % From structure
        xY_sess     = rois{sess}.xY;
        Y_sess      = rois{sess}.Y;
    end
        
    options.avg_method = 'mean';
    
    % Choose voxel-wise or eigenvariate data
    if options.voxel_wise
        y = xY_sess.y;
    else
        switch options.avg_method
            case 'mean'
                y = mean(xY_sess.y, 2);
            case 'eigen'
                y = Y_sess;
            otherwise
                error('Unknown timeseries summary method');
        end
    end
       
    % Validate voxel locations match previous session
    if options.avg_sess
        if sess > 1 && ~all(XYZmm(:) == xY_sess.XYZmm(:))
            error('Could not average: different voxels per session');
        end
        XYZmm = xY_sess.XYZmm;
    end
    
    if ~isfield(Y,'y') || isempty(Y.y)
        Y.y = y;
    elseif options.avg_sess
        % Sum over sessions
        Y.y = Y.y + y;
    else
        % Concatenate
        Y.y = [Y.y; y];
    end
end

if options.avg_sess
    % Complete averaging over sessions
    Y.y = Y.y ./ n_sess;    
end

% check scaling of Y (enforcing a maximum effect size of 4)
scale   = max(max((Y.y))) - min(min((Y.y)));
scale   = 4/max(scale,4);
Y.y     = Y.y*scale;
Y.scale = scale;

% Prepare data
TR   = SPM.xY(1).RT;
Y.dt = TR;
Y.X0 = ones(size(Y.y,1),1);

% Choose microtime offset (secs)
if ~isfield(options,'delays') || isempty(options.delays)
    options.delays = TR / 2;
end

% Convert input timing from seconds -> microtime bins
bins_per_second = 1 / U(1).dt;
bins_per_TR     = bins_per_second * TR;
nscans          = size(Y.y,1);

for t = 1:length(U)
    start_bin = ceil( U(t).ons / U(t).dt) + 1;
    end_bin   = start_bin + (U(t).dur * bins_per_second) - 1;

    U(t).ind   = start_bin : end_bin;           % Microtime bins (from 1)
    U(t).nbins = nscans*bins_per_TR;            % Total bins
end     

% Prepare model
M = specify_model(U,size(Y.y,1),size(Y.y,2),options);

% Package
PRF.xY      = xY_sess;
PRF.M       = M;
PRF.U       = U;
PRF.Y       = Y;
PRF.options = options;

% Get model name if available
if isfield(options,'name')
    PRF.name = ['PRF_' options.name];
elseif isfield(xY_sess,'name') && ~isempty(xY_sess.name)
    PRF.name = ['PRF_' xY_sess.name];
end

save_prf(PRF,SPM.swd);

% -------------------------------------------------------------------------
function PRF = estimate(PRF, est_options)
% Estimate a PRF model
%
% PRF         - model
% est_options - estimation options

% Set default options
try est_options.init;       catch, est_options.init       = 'GLM_P'; end
try est_options.use_parfor; catch, est_options.use_parfor = false; end
try est_options.nograph;    catch, est_options.nograph    = false; end

if ischar(PRF)
    if size(PRF,1) > 1
        error('Please provide only one PRF .mat file for estimation');
    end
    out_dir = fileparts(PRF);
    PRF     = load(PRF);
    PRF     = PRF.PRF;
    do_save = true;
else
    out_dir = pwd;
    do_save = false;
end

ny = size(PRF.Y.y,2);

% Unpack PRF
M     = PRF.M;
U     = PRF.U;
Y     = PRF.Y;
IS    = M.IS;
pE    = M.pE;
pC    = M.pC;

% Results structures
Ep = cell(1,ny);
Cp = cell(1,ny);
Eh = nan(1,ny);
F  = nan(1,ny);

% Voxels to estimate
try
    voxels = est_options.voxels;
catch
    voxels = 1:ny;
end

M.nograph = est_options.nograph;

P = {};

tic
if est_options.use_parfor
    % Run with parallel toolbox
    parfor i = voxels
        if ny > 1, fprintf('Voxel %d of %d\n', i, ny); end
             
        % Initialize priors
        [pE{i},pC{i},P{i}] = initialize_model(M,U,Y.y,i,est_options);
        
        % Model updated with initialized priors
        M2 = M;
        M2.pE = pE{i};
        M2.pC = pC{i};
        M2.P  = P{i};
        
        % Fit        
        [Ep{i},Cp{i},Eh(i),F(i)] = fit_model(M2,U,Y,i);

    end    
else
    % Run single threaded
    for i = voxels 
        if ny > 1, fprintf('Voxel %d of %d\n', i, ny); end

        % Initialize priors
        [pE{i},pC{i},P{i}] = initialize_model(M,U,Y.y,i,est_options);

        % Model updated with initialized priors
        M2 = M;
        M2.pE = pE{i};
        M2.pC = pC{i};
        M2.P  = P{i};
                
        % Fit
        [Ep{i},Cp{i},Eh(i),F(i)] = fit_model(M2,U,Y,i);                           
    end
end
est_time = toc;

% Store initialized priors / starting value
M.pE = pE;
M.pC = pC;
M.P  = P;

% Calculate posterior probability for each parameter
Pp    = cell(1,ny);
for i = 1:ny
    if ~isempty(Ep{i})
        Pp{i} = 1-spm_Ncdf(abs(spm_vec(M.pE{i})),abs(spm_vec(Ep{i})),diag(Cp{i}));
        Pp{i} = spm_unvec(Pp{i},Ep{i});
    end
end

% Pack
PRF.Ep = Ep;
PRF.Cp = Cp;
PRF.Pp = Pp;
PRF.Eh = Eh;
PRF.F  = F;
PRF.Y  = Y;
PRF.M  = M;
PRF.options.estimation = est_options;
PRF.est_time = est_time;

% Add in predicted timeseries (if there's only one timeseries)
if ny == 1 && ~isempty(Ep{1})
    PRF.y = feval(IS,Ep{1},M,U);
end 

% Save
if do_save
    save_prf(PRF,out_dir);
end

% -------------------------------------------------------------------------
function [pE,pC,P] = initialize_model(M,U,y,i,est_options)
% Initializes the parameters of a pRF model based on a rapid search
%
% M - pRF model specification
% U - inputs
% y - data
% i - Voxel index
% est_options - User provided estimation options
%
% Returns:
% pE - Initialized prior expectation of parameters
% pC - Initialized prior covariance of parameters
% P  - Starting values of parameter search

pE = M.pE{i};
pC = M.pC{i};
P  = pE;

switch upper(est_options.init)
    case 'GLM'
        fprintf('Initializing priors using GLM\n');
        
        P_glm = feval(M.IS,pE,M,U,'glm_initialize',y(:,i));
        
        fields = fieldnames(P_glm);
        for i = 1:length(fields)
            pE.(fields{i}) = P_glm.(fields{i});
        end
        
    case 'GLM_P'
        fprintf('Initializing starting value using GLM\n');
        
        P_glm = feval(M.IS,pE,M,U,'glm_initialize',y(:,i));
        
        fields = fieldnames(P_glm);
        for i = 1:length(fields)
            P.(fields{i}) = P_glm.(fields{i});
        end
        
    case 'NONE'
        % No initialization
    otherwise
        error('Unknown initialization option for PRF estimation');
end

% -------------------------------------------------------------------------
function [Ep,Cp,Eh,F] = fit_model(M,U,Y,idx)
% Fit a PRF model to a single timseries (can be run within a parfor loop)
%
% M             - model
% U             - inputs
% XYZmm         - coordinates of each timeseries
% model_options - options from PRF structure
% est_options   - estimation options
%
% Ep            - Posterior means
% Cp            - Posterior variances
% Eh            - Posterior noise log precision
% F             - Free energy

Y.y  = Y.y(:,idx);

if ~isfield(Y,'X0'), Y.X0 = ones(length(Y.y),1); end

% Fit model
try
    [Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);
catch e
    warning('PRF failed to converge');
    Ep = [];
    Cp = [];
    Eh = NaN;
    F  = NaN;
end

% -------------------------------------------------------------------------
function M = specify_model(U,ns,nv,options)
% Specify a PRF model
%
% U       - inputs
% ns      - number of scans
% nv      - number of voxels
% options - model specification options
%
% M       - model

% Response function
try options.model; catch, options.model = 'spm_prf_fcn_gaussian'; end
M.IS = options.model;

% Parameter scaling
try 
    M.pmax = U(1).pmax;
    M.pmin = U(1).pmin;
catch
    M.pmax = 2;
    M.pmin = 0;
end

% Neuronal priors
[pE,pC] = feval(M.IS,[],M,[],'get_priors');

% Neurovascular priors
pv = 1/128;
pE.transit = 0;
pE.decay   = 0;
pC.transit = pv;
pC.decay   = pv;    

% Field strength
M.B0 = options.B0;

% BOLD prior - Field-dependent fixed epsilon (Heinzle et al. 2016)
switch M.B0
    case 1.5
        pE.epsilon = 0.25;
        pC.epsilon = 0.04;
    case 3
        pE.epsilon = -0.78;
        pC.epsilon = 0.24;
    case 7
        pE.epsilon = -3.99;
        pC.epsilon = 0.83;
    otherwise
        error('Please choose your closest field strength - 1.5,3 or 7');
end

% Apply to all PRFs
for v = 1:nv
    M.pE{v} = pE;
    M.pC{v} = pC;
end

% Starting values
if isfield(options,'P') && ~isempty(options.P)
    M.P = options.P;
end

% Model spec
M.f  = 'spm_prf_fx';
M.g  = 'spm_prf_gx';
M.x = zeros(4,1);   % Initial hemodynamic state
M.m = 1;            % Number of inputs
M.n = length(M.x);  % Number of states
M.l = 1;            % Number of outputs
M.N = 32;
M.dt = U(1).dt;
M.ns = ns;
M.TE = options.TE;
M.hE = options.hE;
M.T0 = max(ceil(options.delays / M.dt),1);

% -------------------------------------------------------------------------
function save_prf(PRFs, parent_folder)
% Save one or more PRFs
%
% PRF           - A single PRF or a [1 x n] cell array of PRFs
% parent_folder - Where to save the model

if ~iscell(PRFs), PRFs = {PRFs}; end
if nargin < 2, parent_folder = pwd; end

if ~exist(parent_folder,'file') && ~isempty(parent_folder)
    mkdir(parent_folder);
end

for r = 1:length(PRFs)
    PRF = PRFs{r};
    
    % Get PRF name    
    if isfield(PRF,'name') && ~isempty(PRF.name)
        name = PRF.name;
    else
        name = ['PRF_R' num2str(r)];
    end
    
    % Save
    PRF.name = name;    
    filename = fullfile(parent_folder, [name '.mat']);    
    save(filename,'PRF');
end

% -------------------------------------------------------------------------
function split_model(PRF,num_models,out_dir)
% Split a PRF file containing multiple timeseries into several files
%
% PRF        - original model
% num_models - target number of model files
% out_dit    - output directory

PRF0 = PRF;

nv = size(PRF0.Y.y,2);

voxels_per_model = ceil(nv / num_models);

first_voxel       = 1:voxels_per_model:nv;
last_voxel        = first_voxel + voxels_per_model - 1;
last_voxel(end)   = nv;
num_models        = length(first_voxel);

PRF0.xY = rmfield(PRF0.xY,'y');

for p = 1:num_models
    fprintf('Writing model %d\n',p);
    vox_idx = first_voxel(p):last_voxel(p);
    
    PRF = extract_model(PRF0,vox_idx);
    
    PRF.name = sprintf('%s_%d', PRF0.name, p);
    save_prf(PRF, out_dir);
end
% -------------------------------------------------------------------------
function merge_models(P,out_dir)
% Merge multiple PRF model files into one
%
% P       - cell array of filename strings
% out_dir - output directory

XYZmm = [];
Y.y = [];

pE = {};
pC = {};
Ep = {};
Cp = {};
Pp = {};
Eh = [];
F  = [];
est_time = [];
nm = []; % Number of PRF models in each file

for m = 1:length(P)
    fprintf('Reading PRF file %d\n',m);
    PRF=load(P{m});
    PRF=PRF.PRF;
    
    XYZmm = [XYZmm PRF.xY.XYZmm];
    Y.y   = [Y.y PRF.Y.y];
    
    pE = [pE PRF.M.pE];
    pC = [pC PRF.M.pC];
    Ep = [Ep PRF.Ep];
    Cp = [Cp PRF.Cp];
    Pp = [Pp PRF.Pp];
    Eh = [Eh PRF.Eh];
    F  = [F  PRF.F];   
    
    est_time = [est_time PRF.est_time];
    
    nm(m) = length(PRF.Ep);
end

name = PRF.name;

% Remove trailing underscore
idx = strfind(PRF.name,'_');
if ~isempty(idx)
    idx  = idx(end);
    name = name(1:idx-1);
end

name = [name '_combined'];

% Pack & save PRF
PRF.xY.XYZmm = XYZmm;
PRF.Y.y      = Y.y;
PRF.M.pE     = pE;
PRF.M.pC     = pC;
PRF.Ep       = Ep;
PRF.Cp       = Cp;
PRF.Pp       = Pp;
PRF.Eh       = Eh;
PRF.F        = F;
PRF.name     = name;
PRF.est_time = est_time ./ nm;

save_prf(PRF,out_dir);
% -------------------------------------------------------------------------
function PRF = extract_model(PRF0,vox_idx)
% Extract a single PRF from a multi-PRF file

PRF = PRF0;
if isfield(PRF,'Ep')
    PRF.Ep = PRF0.Ep(vox_idx);
    PRF.Cp = PRF0.Cp(vox_idx);
    PRF.Pp = PRF0.Pp(vox_idx);
    PRF.Eh = PRF0.Eh(vox_idx);
    PRF.F  = PRF0.F(vox_idx);
end

if isfield(PRF,'M') && isfield(PRF.M,'pE')
    PRF.M.pE = PRF0.M.pE(vox_idx);
    PRF.M.pC = PRF0.M.pC(vox_idx);    
end

PRF.Y.y = PRF0.Y.y(:,vox_idx);

if isfield(PRF0.xY,'XYZmm')
    PRF.xY.XYZmm = PRF0.xY.XYZmm(:,vox_idx);
end
% -------------------------------------------------------------------------
function Ep = get_corrected_parameters(PRF)
% Returns parameters corrected by converting from latent variables

Ep = cell(1,length(PRF.Ep));
for v = 1:length(PRF.Ep)
    Ep{v} = feval(PRF.M.IS, PRF.Ep{v}, PRF.M, PRF.U, 'get_parameters');
end