function varargout = spm_prf_review(varargin)
% Review a PRF model. If only one timeseries is included, parameter plots
% are shown. Otherwise, an 3D parameter map is displayed. If the parameter
% map doesn't yet exist, it will be generated.
%
% -------------------------------------------------------------------------
% Review the given model:
%
% spm_prf_review(PRF,idx)
%
% PRF - model
% idx - (Optional). Voxel index within PRF. If blank, generates and 
%       displays a results image.
%
% -------------------------------------------------------------------------
% Generate results image only:
%
% spm_prf_review(PRF,'generate_images',SPM)
%
% PRF - model
% SPM - SPM structure
%
% -------------------------------------------------------------------------
% View a posterior probability map
%
% spm_prf_review(PRF,'view_comparison',filename)
%
% PRF - PRF
% filename - PPM image filename
%
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

% Validate inputs
try 
    PRF = varargin{1};
catch
    [PRF,sts] = spm_select(1,'mat','Please select PRF file');
    if ~sts, return; end
end

% Record source directory
if ischar(PRF)
    prf_dir = fileparts(PRF);        
else
    prf_dir = pwd;
end

% Load PRF
if ischar(PRF)
    PRF = load(PRF);
    PRF = PRF.PRF;   
end

% Store source directory in PRF
if ~isfield(PRF,'dir') || isempty(PRF.dir)
    PRF.dir = prf_dir;
end

% Check estimated
if ~isfield(PRF,'Ep')
    error('Please estimate the model before reviewing');
end

if length(varargin) >= 2 && ischar(varargin{2})
    % An action string has been provided
    action = varargin{2};
    
    switch action
        case 'generate_images'
            % Generate results image only

            % Load SPM
            SPM = load_spm(PRF);
            if isempty(SPM); return; end

            % Generate
            varargout{1} = ...
                fullfile(PRF.dir,create_results_img(PRF,SPM.swd,SPM.VM,true));
        case 'view_comparison'
            % Get results image
            fname = varargin{3};
                                   
            % Create figure
            spm_figure('GetWin','Graphics');
            spm_clf;    
            clear global st                           
            
            % Display
            review_prf_map(PRF,[],fname,false);
        otherwise
            error('Unknown action');
    end
else
    % View PRF results
    
    % Check for PRF voxel index in second argument
    if length(varargin) >= 2 && ~isempty(varargin{2})
        idx = varargin{2};
    else
        idx = [];
    end
    
    % Check for mask in third argument
    if length(varargin) >= 3
        xY = varargin{3};
        if ischar(xY)
            xY = load(xY);
            xY = xY.xY;
        end
    else
        xY = [];
    end
    
    % Create figure
    spm_figure('GetWin','Graphics');
    spm_clf;    
    clear global st

    if length(PRF.Ep) == 1 || ~isempty(idx)
        % Report results for a single timeseries  
        if isempty(idx) 
            idx = 1;
        end
        idx = max(idx,1);
        review_single_prf(PRF, idx);        
    else
        % Report results across an image
        review_prf_map(PRF,xY);      
        
        % Diagnostic text
        n_failed = sum(cellfun(@isempty,PRF.Ep));
        fprintf('Timeseries with failed estimation: %d\n', n_failed);

    end
    
end
% -------------------------------------------------------------------------
function review_prf_map(PRF,xY,fname,is_model_parameter)
% Displays a brain map of PRF parameters or model comparison results
%
% PRF                - the PRF structure
% xY                 - mask
% fname              - Filename of the image to display. If not provided,  
%                      it will be generated.
% is_model_parameter - if true, the data being displayed are model
%                      parameters. If false, they are model comparison 
%                      results.

% Load SPM
SPM = load_spm(PRF);
if isempty(SPM); return; end                      

% Get or create xPRF structure
try
    xPRF = evalin('base','xPRF');
catch                               
    xPRF = struct('PRF',PRF,...
                  'voxel_idx',0,...
                  'param_idx',1,...
                  'selected_voxel_idx',-1,...
                  'V',SPM.Vbeta(1),...
                  'xY',xY);
end

xPRF.xY = xY;

if nargin < 3 || is_model_parameter
    % Overwrite if a mask is given
    force_overwrite = ~isempty(xY);

    if isfield(PRF,'pimg') && (exist(fullfile(PRF.dir,PRF.pimg),'file') > 0) ...
            && ~force_overwrite && isempty(xPRF.xY)
        % Get existing parameter map
        fname = PRF.pimg;
    else
        % Create parameter map
        fname    = create_results_img(PRF,SPM.swd,xPRF.V,true,xY);
        PRF.pimg = fname;
        
        % Save PRF with map filename for next time
        if isempty(xPRF.xY)            
            if isfield(PRF,'dir')
                save(fullfile(PRF.dir,[PRF.name '.mat']),'PRF');
            end
        end
    end
    
    % Prepend path
    fname = fullfile(PRF.dir,fname);
    
    % Get index in 4D parameter map nifti
    fname = [fname ',' num2str(xPRF.param_idx)];        
end

% Get indices of non-nan voxels
V = spm_vol(fname);
Y = spm_read_vols(V);
is_nz_vox = ~isnan(Y(:));

if nargin < 4
    is_model_parameter = true;
end

% Background image
if isfield(xPRF,'anatomical')
    bg = xPRF.anatomical;
else
    bg = fullfile(SPM.swd,SPM.VM.fname);
end

% Add drop-down menu
spm_figure('GetWin','Graphics');
h = uicontrol('Style','Popupmenu','Units','normalized', 'Position',[0.1 0.89 0.5 0.1],...
    'Tag','spm_image:window');

% Add parameter selector
if is_model_parameter
   fields = spm_fieldindices(PRF.Ep{1},1:length(spm_vec(PRF.Ep{1})));
   
   set(h,'Callback',@selected_param_changed, ...
         'ToolTipString','Parameter to display',...
         'String',fields,...
         'Value',xPRF.param_idx);
else
   set(h,'Callback','', ...
         'ToolTipString','',...
         'String','Posterior probability map',...
         'Value',1,'Enable','off');
end

% Add anatomical button
uicontrol('Style','PushButton','Units','normalized', ...
          'Position',[0.62 0.96 0.15 0.03],'String','Anatomical',...
          'Callback',@choose_anatomical);
      
% Add surface button
uicontrol('Style','PushButton','Units','normalized', ...
          'Position',[0.8 0.96 0.15 0.03],'String','Surface',...
          'Callback',@project_surface);

% Colormap
c = jet();

% Orthogonal views
spm_orthviews('Image',bg,[0 0.59 1 0.4]);
spm_orthviews('Addtruecolourimage',1,V,c,0.5);
spm_orthviews('AddContext');
spm_orthviews('Redraw');

global st

% Store PRF voxel locations as indices in displayed image
% mm -> vox
mat = V.mat;
dim = V.dim;
XYZmm = PRF.xY.XYZmm;
XYZ = mat \ [XYZmm; ones(1,size(XYZmm,2))];
XYZ = round(XYZ);
% vox -> idx
idx = sub2ind(dim,XYZ(1,:),XYZ(2,:),XYZ(3,:));        

% Axes for plots
ax = uipanel('Units','normalized','Position',[0 0 1 0.55],'BackgroundColor',[1 1 1],'BorderWidth',0);

% Update xPRF
xPRF.ax = ax;
xPRF.idx = idx;
xPRF.PRF = PRF;
xPRF.is_nz_vox = is_nz_vox;
if is_model_parameter
    xPRF.sel_map_image = '';
else
    xPRF.sel_map_image = fname;
end
assignin('base','xPRF',xPRF);

% Add listeners
st.callback = @selected_voxel_changed;
% -------------------------------------------------------------------------
function selected_param_changed(varargin)
% Callback after the user selects a different parameter to view
p   = get(varargin{1},'Value');

xPRF = evalin('base','xPRF');
xPRF.param_idx = p;
assignin('base','xPRF',xPRF);

% spm_prf_review(xPRF.PRF,...
%                xPRF.voxel_idx,...
%                xPRF.param_idx,...
%                xPRF.selected_voxel_idx);
spm_prf_review(xPRF.PRF);
% -------------------------------------------------------------------------
function selected_voxel_changed(varargin)
% Callback after the user moves the crosshairs in the display
global st

xPRF = evalin('base','xPRF');

mat = xPRF.V.mat;
dim = xPRF.V.dim;

%mat    = st.vols{1}.mat;
%dim    = st.vols{1}.dim;
centre = st.centre;

parent = xPRF.ax;

% Selected mm -> vox
XYZ = mat \ [centre; 1];
XYZ = round(XYZ);

% vox -> idx
selected_idx = sub2ind(dim,XYZ(1),XYZ(2),XYZ(3));

% Record selected voxel
idx = find(xPRF.idx == selected_idx);
fprintf('Selected PRF %d\n', idx);
xPRF.selected_voxel_idx = idx;
assignin('base','xPRF',xPRF);

% Clear display
h = allchild(xPRF.ax);
delete(h);

% Check if that voxel was significant and in-mask, if so display
if ~isempty(idx) && xPRF.is_nz_vox(selected_idx)
    review_single_prf(xPRF.PRF,idx,parent);
end

% -------------------------------------------------------------------------
function choose_anatomical(varargin)

[t,sts] = spm_select(1,'image','Please select anatomical');
if ~sts, return; end

% Update
xPRF = evalin('base','xPRF');
xPRF.anatomical = t;
assignin('base','xPRF',xPRF);

% Refresh
spm_prf_review(xPRF.PRF);

% -------------------------------------------------------------------------
function project_surface(varargin)
% Callback for projection to the surface button

xPRF = evalin('base','xPRF');
PRF  = xPRF.PRF;
pidx = xPRF.param_idx;

% Map image (e.g. model comparison results) if selected
sel_map_image = xPRF.sel_map_image;

spm_prf_review_surface(PRF,pidx,sel_map_image);

% -------------------------------------------------------------------------
function review_single_prf(PRF,idx,parent)
% Creates review window with plots for a single timeseries' results
%
% PRF    - model
% idx    - voxel index within PRF. If blank, displays results image
% parent - parent UI element

spm_figure('GetWin','Graphics');

if nargin < 3 || isempty(parent)
    parent = gcf;
end

% Plot settings
rows = 3;
cols = 4;
titlesize = 16;
legendsize = 12;

% List X,Y stimulus coordinates
pmax = PRF.M.pmax/2;
x_bins = -pmax:0.5:pmax;
y_bins = x_bins;
b      = length(x_bins);
[x2,y2] = meshgrid(x_bins,y_bins);
xy = [x2(:) y2(:)];

% PRF prior response
z_prior  = feval(PRF.M.IS, PRF.M.pE{idx}, PRF.M, PRF.U, 'get_response', xy);
z_prior  = reshape(z_prior,b,b);

% PRF with fixed beta parameters for sampling
PRF_tmp = PRF;
fields  = fieldnames(PRF_tmp.M.pE{idx});
for f = 1:length(fields)
    if strcmp(fields{f}(1:4),'beta')
        PRF_tmp.M.pC{idx}.(fields{f}) = 0;
    end
end

% Prior and posterior predictive densities (fixed prior beta parameter)
[x,x,prior_pd,post_pd] = spm_prf_get_ppd(PRF_tmp,xy,idx,1000);
prior_pd = reshape(prior_pd,b,b);
post_pd  = reshape(post_pd,b,b);
clear PRF_tmp;

% Priors
pE = PRF.M.pE{idx};
pC = PRF.M.pC{idx};

is_estimated = ~isempty(PRF.Ep{idx});
if is_estimated
    % PRF posterior response
    z_post   = feval(PRF.M.IS, PRF.Ep{idx}, PRF.M, PRF.U, 'get_response', xy);
    z_post   = reshape(z_post,b,b);
    
    % Get parameters
    Ep = PRF.Ep{idx};
    Cp = PRF.Cp{idx};

    % Get predicted timeseries if not saved
    if ~isfield(PRF,'y') || isempty(PRF.y)
        [y,Z] = feval(PRF.M.IS,PRF.Ep{idx},PRF.M,PRF.U);
        PRF.y = y;   % BOLD
        PRF.x = Z.u; % Neuronal
    end
    
    % Calculate residuals
    if ~isfield(PRF,'R')
        Y = PRF.Y;
        R = Y.y(:,idx) - PRF.y;
        PRF.R = R - Y.X0*spm_inv(Y.X0'*Y.X0)*(Y.X0'*R);
    end
end

% Sample to get underlying parameters for bar plot
if true
    
    % Sample prior and posterior (unlike above, we don't fix beta)
    [spE,sEp] = spm_prf_get_ppd(PRF,xy,idx,1000);

    % Correct sampled priors
    for s = 1:size(spE,2)        
        P        = spm_unvec(spE(:,s),PRF.M.pE{idx});
        spE(:,s) = spm_vec(feval(PRF.M.IS, P, PRF.M, PRF.U, 'get_parameters'));         
    end
    
    % Correct sampled posteriors
    if is_estimated        
        for s = 1:size(sEp,2)        
            P        = spm_unvec(sEp(:,s),PRF.Ep{idx});
            sEp(:,s) = spm_vec(feval(PRF.M.IS, P, PRF.M, PRF.U, 'get_parameters'));         
        end
    end    
end

% Calculate explained variance
PSS   = sum(PRF.y.^2);
RSS   = sum(PRF.R.^2);
exp_var = 100*PSS/(PSS + RSS);

% Plot PRFs
if is_estimated, num_plots = 4; else num_plots = 2; end

for i = 1:num_plots
    if i > 2 && ~is_estimated
        break;
    end
    
    subplot(rows,cols,i,'Parent',parent);
    colormap(jet);
    switch i
        case 1
            imagesc(z_prior);
            title('Prior','FontSize',titlesize);
        case 2 
            imagesc(prior_pd);
            title('Prior PD','FontSize',titlesize);
        case 3
            imagesc(z_post);
            title('Posterior','FontSize',titlesize);
        case 4
            imagesc(post_pd);
            title('Post PD','FontSize',titlesize);           
    end
    
    set(gca,'YDir','normal','XTick',[1 ceil(b/2) b],'YTick',[1 ceil(b/2) b]);
    set(gca,'XTickLabel',[-pmax 0 pmax],'YTickLabel',[-pmax 0 pmax]); 
    set(gca,'FontSize',legendsize); 
    axis square;    
end

if ~is_estimated, return; end;
 
% Median plot
subplot(rows,cols,5:8);

jitter     = 0.09;
x          = 1:size(spE,1);
fields     = spm_fieldindices(PRF.Ep{idx},1:length(spm_vec(PRF.Ep{idx})));
hpE        = spm_plot_quartiles(spE',x-jitter,'Color',[0.6 0.6 0.6]); hold on;
[hEp,hBar] = spm_plot_quartiles(sEp',x+jitter,'Color',[0 0 0]);

set(gca,'XTickLabel',fields,'XTick',1:length(fields),'FontSize',10);
title('Parameters','FontSize',titlesize);
legend([hpE(1) hEp(1) hBar(1)],{'Prior','Posterior','IQR'},'FontSize',12);

% Estimated timeseries
subplot(rows,cols,9:12);
x_response = (1:size(PRF.Y.y,1))*PRF.Y.dt;
plot(x_response,PRF.y+PRF.R,'Color',[0.3 0.3 0.3],'LineStyle','--'); hold on;
plot(x_response,PRF.y,'Color','r');
hold off;
xlabel('Time (secs)','FontSize',legendsize);
legend({'BOLD','Model'},'Location','southeast');
title(sprintf('Outputs - exp. var = %1.0f%%',exp_var),'FontSize',titlesize);
set(gca,'FontSize',legendsize);


% F
pC = PRF.M.pC{idx};
pE = PRF.M.pE{idx};
Cp = PRF.Cp{idx};
Ep = PRF.Ep{idx};

rE = spm_vec(pE);
rC = spm_vec(pC);

% Indices of parameters to keep TODO move this
np = length(rE);
q  = zeros(1,np);
q(end-3:end) = 1;

rC(q~=1)=0;

[F,sE,sC] = spm_log_evidence(Ep,Cp,pE,diag(spm_vec(pC)),rE,diag(rC));
disp(F);

% Convert to post. prob - from spm_dcm_compare
F    = [F 0];
i    = F < (max(F) - 32);
P    = F;
P(i) = max(F) - 32;
P    = P - min(P);
P    = exp(P);
P    = P/sum(P);

fprintf('Reduced (Null) - Full model. log BF: %2.2f, P(null): %2.2f\n',F(1),P(1));

% -------------------------------------------------------------------------
function out_file = create_results_img(PRF,out_dir,V,overwrite,xY)
% Create results images, one per parameter, from a PRF
%
% PRF       - estimated model
% 
% overwrite - tf whether to overwrite existing file
%
% out_file - created 4D nifti

% Filename for parameter file
out_file = [PRF.name '.nii'];
fname    = fullfile(out_dir,out_file);

if exist(fname,'file') && ~overwrite
    return;
end

if ~isfield(PRF,'Ep')
    error('Please estimate this model before continuing');
end

disp('Creating results image');

nv = length(PRF.Ep);
np = length(spm_vec(PRF.Ep{1}));

% Identify significant voxels to display
nv  = length(PRF.Ep);
sig_idx = zeros(1,nv);
for v = 1:nv
    if isempty(PRF.Ep{v})
        warning('Voxel %d not estimated', v);
        continue;
    end
    sig_idx(v) = feval(PRF.M.IS, PRF.Ep, PRF.M, PRF.U, 'is_above_threshold', PRF.Cp, v, 0.95);
end
sig_idx = find(sig_idx);

fprintf('%d out of %d voxels were above threshold\n',length(sig_idx),nv);

% Spatial mask
if nargin >=5 && ~isempty(xY)
    [~, ~, included_voxels] = spm_ROI(xY, PRF.xY.XYZmm);
    is_included = ismember(sig_idx, included_voxels);
    sig_idx     = sig_idx(is_included == 1);
end

% Collate parameters and confidence intervals
Ep_all = nan(nv,np);
Ci_all = nan(nv,np);
ci     = spm_invNcdf(1 - 0.05);
for v = sig_idx
    % Correct parameters for display
    Ep = feval(PRF.M.IS, PRF.Ep{v}, PRF.M, PRF.U, 'get_parameters');

    % Get 95% CI
    Ci_all(v,:) = ci * sqrt(diag(PRF.Cp{v}));
        
    % Store
    Ep_all(v,:) = spm_vec(Ep);
end

Ep_all(Ep_all == 0) = nan;

% mm -> vox
Q   = ones(1,size(PRF.xY.XYZmm,2));
XYZ = V.mat \ [PRF.xY.XYZmm; Q];
XYZ = XYZ(1:3,:);

% Deal with rounding errors - constrain between 1 and image dimension
XYZ(1,:) = min(XYZ(1,:), V.dim(1));
XYZ(2,:) = min(XYZ(2,:), V.dim(2));
XYZ(3,:) = min(XYZ(3,:), V.dim(3));
XYZ(1,:) = max(XYZ(1,:), 1);
XYZ(2,:) = max(XYZ(2,:), 1);
XYZ(3,:) = max(XYZ(3,:), 1);

% vox -> idx
j = sub2ind(V.dim,XYZ(1,:),XYZ(2,:),XYZ(3,:));
j = round(j);

start_dir = pwd;
cd(out_dir);

% Make parameter images
fprintf('Writing parameter maps\n');
for p = 1:np
    Y = nan(V.dim);
    Y(j) = Ep_all(:,p);
    
    V.n(1)  = p;
    V.fname = out_file;
    spm_write_vol(V,Y);
end

% Make confidence interval images
fprintf('Writing confidence interval maps\n');
for p = 1:np
    % Upper
    Y = nan(V.dim);
    Y(j) = Ep_all(:,p) + Ci_all(:,p);
    
    V.n(1)  = p;
    V.fname = ['CI_upper_' out_file];
    spm_write_vol(V,Y);
    
    % Lower
    Y = nan(V.dim);
    Y(j) = Ep_all(:,p) - Ci_all(:,p);
    
    V.n(1)  = p;
    V.fname = ['CI_lower_' out_file];
    spm_write_vol(V,Y);    
end

cd(start_dir);

% -------------------------------------------------------------------------
function SPM = load_spm(PRF)
% Load SPM from the same folder as the PRF, or prompt for an SPM
%
% PRF - PRF structure
% SPM - SPM structure

SPM = [];
try
    SPM = load(fullfile(PRF.dir,'SPM.mat'));        
    spm_dir = PRF.dir;
catch
    [P,sts] = spm_select(1,'mat','Please select SPM.mat',{},pwd,'SPM.mat');
    if sts
        SPM = load(P);
        spm_dir = fileparts(P);
    else
        return;
    end
end    

SPM = SPM.SPM;

if isempty(SPM.swd)
    SPM.swd = spm_dir;
end