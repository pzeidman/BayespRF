function spm_prf_editor(PRF,mask,posterior)
% Tool for exploring the prior and posterior settings of a pRF model
%
% PRF       - the pRF model, which may or may not be estimated
% mask      - if set, this is a 2D binary stimulus matrix which is overlaid
%             on the pRF response
% posterior - if true, the posterior parameters are used rather than the
%             priors [default: false]
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

% Load PRF
if nargin < 1
    [P,sts] = spm_select(1,'mat','Select PRF',{},pwd,'^PRF_.*mat$');
    if ~sts, return; end
    PRF = load(P);
    PRF = PRF.PRF;
end

% Just use first voxel
if isfield(PRF,'Ep') && iscell(PRF.Ep)
    PRF.Ep = PRF.Ep{1};
    PRF.Cp = PRF.Cp{1};
end
if isfield(PRF.M,'pE') && iscell(PRF.M.pE)
    PRF.M.pE = PRF.M.pE{1};
    PRF.M.pC = PRF.M.pC{1};
end

if nargin < 2
    mask = [];
end

% Default to prior if not estimated
if ~isfield(PRF,'Ep')
    posterior = false;
end

% Chose posterior or prior params
if nargin >= 3 && posterior
    Ep = PRF.Ep;
    Cp = PRF.Cp;
else
    Ep = PRF.M.pE;
    Cp = PRF.M.pC;
end

% Check covariance is a matrix
if isstruct(Cp)
    Cp = diag(spm_vec(Cp));
end

% Create figure
f = figure('Name','PRF Designer','NumberTitle','Off');
p = get(f,'Position');
p(2) = 50;
p(3) = 700;
p(4) = 700;
set(f,'Position',p);

% Add controls
handles = add_controls(f);

% Set default display options
options = struct('nsamp',1000,'disp_ppd',true);

% Store in figure
data = struct();
data.handles = handles;
data.Ep      = Ep;
data.Cp      = Cp;
data.PRF     = PRF;
data.mask    = mask;
data.options = options;
data.sel_param_idx = 1;
set(gcf,'UserData',data);

% Populate controls
update_controls(Ep,Cp,options,handles);

display_PRF;

% -------------------------------------------------------------------------
function update_controls(Ep,Cp,options,h)
% Update controls to given values
%
% Ep - Structure of expected values
% Cp - Covariance matrix
% options - structure with disp_ppd and nsamp
% h  - Handles structure from add_controls()

% Format for displaying floating point numbers
frmt = '%3.3f';

% Get structure array of variances
Vp = diag(Cp);
Vp = spm_unvec(Vp,Ep); %#ok<NASGU>

% Display parameter names
pnames = spm_fieldindices(Ep,1:length(spm_vec(Ep)));
set(h.pselect,'String',pnames);

% Identify selected parameter
selected_idx   = get(h.pselect,'Value');
selected_pname = pnames{selected_idx};

% Display mean and variance
Ep_p = eval(['Ep.' selected_pname]);
Vp_p = eval(['Vp.' selected_pname]);
set(h.exp_edit,'String',sprintf(frmt,Ep_p));
set(h.var_edit,'String',sprintf(frmt,Vp_p));

% Set sliders
set(h.exp_slider,'Min',Ep_p - 3,'Max',Ep_p + 3,'Value',Ep_p);
set(h.var_slider,'Min',0,'Max',Vp_p + 3,'Value',Vp_p);

% Display options
set(h.ppd_samp,'String',options.nsamp);
set(h.ppd_tog,'Value',options.disp_ppd);

% -------------------------------------------------------------------------
function handles = add_controls(f)

% Get figure size
p = get(f,'Position');
w = p(3); h = p(4);

% Add left panel
panel = uipanel('Parent',f,'Units','pixels','Position',[0 0 (w/2) h]);

% Add right axes
m     = 5; % margin
prf_h = (h/3) - (m*3);
prf_w = prf_h;
prf_l = (w/2) + (prf_w/4);

bold_axes = axes('Parent',f,'Units','pixels','Position',[prf_l (prf_h*2)+(m*10) prf_w*1.2 prf_h*0.8],'Color','w');
prf_axes  = axes('Parent',f,'Units','pixels','Position',[prf_l prf_h+m prf_w prf_h],'Color','w');
ppd_axes  = axes('Parent',f,'Units','pixels','Position',[prf_l m       prf_w prf_h],'Color','w');

% Get panel size
p = get(panel,'Position');
w = p(3); h = p(4);

% Add parameter selector
margin     = 5;
listheight = 70;
listbottom = (h-listheight-margin);
pselect = uicontrol('Style','listbox','Parent',panel,'Units','pixels',...
    'Position',[margin listbottom w-(margin*2) listheight],...
    'FontSize',12,'Callback',@pselect_callback);

% Add expectation panel
exp_panel_h = 120;
exp_panel_w = w-(margin*2);
exp_panel_b = listbottom-(margin*2)-exp_panel_h;
exp_panel = uipanel('Parent',panel,'Units','pixels',...
    'Position',[margin exp_panel_b exp_panel_w exp_panel_h],...
    'Title','Expected value','FontSize',12);

% Add expectation editor and slider
exp_edit_w = 80;
exp_edit_h = 30;
exp_edit_b = exp_panel_h-exp_edit_h-(margin*4);
exp_edit = uicontrol('style','edit','Parent',exp_panel,'Units','pixels',...
    'Position',[margin exp_edit_b exp_edit_w exp_edit_h],'FontSize',12,...
    'Callback',@editor_callback);
exp_slider = uicontrol('style','slider','Parent',exp_panel,'Units','pixels',...
    'Position',[(margin*2)+exp_edit_w exp_edit_b-margin (exp_panel_w-exp_edit_w-(margin*4)) exp_edit_h],...
    'Callback',@slider_callback);

% Add confidence interval labels
ci_text_h = 50;
ci_text = uicontrol('style','text','Parent',exp_panel,...
    'Position',[margin exp_edit_b-ci_text_h-(margin*3) exp_panel_w-(margin*2) ci_text_h],...
    'HorizontalAlignment','left','FontSize',12);

% Add variance panel
var_panel_height = 60;
var_panel_b = exp_panel_b-(margin*2)-var_panel_height;
var_panel = uipanel('Parent',panel,'Units','pixels',...
    'Position',[margin var_panel_b exp_panel_w var_panel_height],...
    'Title','Variance (uncertainty)','FontSize',12);

% Add variance editor and slider
var_edit_b = var_panel_height-exp_edit_h-(margin*4);
var_edit = uicontrol('style','edit','Parent',var_panel,'Units','pixels',...
    'Position',[margin var_edit_b exp_edit_w exp_edit_h],'FontSize',12,...
    'Callback',@editor_callback);
var_slider = uicontrol('style','slider','Parent',var_panel,'Units','pixels',...
    'Position',[(margin*2)+exp_edit_w var_edit_b-margin (exp_panel_w-exp_edit_w-(margin*4)) exp_edit_h],...
    'Callback',@slider_callback);

% Add display settings panel
dis_panel_h = 100;
dis_panel_b = var_panel_b-(margin*2)-dis_panel_h;
dis_panel   = uipanel('Parent',panel,'Units','pixels',...
    'Position',[margin dis_panel_b exp_panel_w dis_panel_h],...
    'Title','Display settings','FontSize',12);

% Add PPD toggle
ppd_tog_w = 120;
ppd_tog_h = 30;
ppd_tog_b = dis_panel_h-ppd_tog_h-(margin*4);
ppd_tog   = uicontrol('style','checkbox','Parent',dis_panel,'Units','pixels',...
    'Position',[margin ppd_tog_b ppd_tog_w ppd_tog_h],'String','Display PD',...
    'FontSize',12,'Callback',@display_PRF);

% Add PPD #samples editor and label
samples_h = 30;
samples_w = 80;
samples_b = ppd_tog_b-samples_h;
samples_edit = uicontrol('style','edit','Parent',dis_panel,'Units','pixels',...
    'Position',[margin samples_b samples_w samples_h],'FontSize',12,...
    'Callback',@display_PRF);
uicontrol('style','text','Parent',dis_panel,'Units','pixels',...
    'Position',[margin*2+samples_w samples_b-5 samples_w samples_h],...
    'String','samples','HorizontalAlignment','left','FontSize',12);

% Add reset button
% reset_w = 60;
% reset_h = 20;
% reset_b = samples_b-reset_h-(margin*3);
% reset_btn = uicontrol('style','pushbutton','Parent',dis_panel,'Units','pixels',...
%     'Position',[margin reset_b reset_w reset_h],'String','Reset all',...
%     'FontSize',12);

% Package up handles
handles.bold_axes  = bold_axes;
handles.prf_axes   = prf_axes;
handles.ppd_axes   = ppd_axes;
handles.pselect    = pselect;
handles.exp_edit   = exp_edit;
handles.exp_slider = exp_slider;
handles.ci_text    = ci_text;
handles.var_edit   = var_edit;
handles.var_slider = var_slider;
handles.ppd_tog    = ppd_tog;
handles.ppd_samp   = samples_edit;
%handles.reset_btn  = reset_btn;

% -------------------------------------------------------------------------
function [Ep,Cp,options] = update_from_gui(f)
% Get PRF parameters updated by GUI

data  = get(f,'UserData');
h     = data.handles;
Ep    = data.Ep;
Cp    = data.Cp;
options = data.options;
selected_idx = data.sel_param_idx;

% Identify selected parameter
pnames = get(h.pselect,'String');
selected_pname = pnames{selected_idx};

% Get selected values
sel_exp = str2double(get(h.exp_edit,'String'));
sel_var = str2double(get(h.var_edit,'String'));
eval(sprintf('Ep.%s=%3.3f', selected_pname, sel_exp));
eval(sprintf('Cp(%d,%d)=%3.3f', selected_idx, selected_idx, sel_var));

% Get options
options.nsamp    = str2double(get(h.ppd_samp,'String'));
options.disp_ppd = get(h.ppd_tog,'Value');

% -------------------------------------------------------------------------
function display_PRF(varargin)

% TODO - get figure from control if set
f = gcf;

% Get mask
data  = get(f,'UserData');
mask  = data.mask;
PRF   = data.PRF;
handles = data.handles;
selected_idx = data.sel_param_idx;

% Get parameter settings from GUI
[Ep,Cp,options] = update_from_gui(f);

% Calculate 95% confidence interval
Ep_p = spm_vec(Ep);
Ep_p = Ep_p(selected_idx);
Vp_p = Cp(selected_idx,selected_idx);
ci   = spm_invNcdf(1 - 0.05);
c    = ci * sqrt(Vp_p);
ci_lower = Ep_p - c;
ci_upper = Ep_p + c;
ci_str = sprintf('95%% confidence interval %3.3f to %3.3f',ci_lower,ci_upper);

% Latent variable -> underlying parameter
Ep_exp = feval(PRF.M.IS, Ep, PRF.M, PRF.U, 'get_parameters');
Ep_exp = spm_vec(Ep_exp);
Ep_exp = Ep_exp(selected_idx);

% Latent variable upper CI -> underlying parameter
Ep_up  = spm_vec(Ep);
Ep_up(selected_idx) = ci_upper;
Ep_up  = spm_unvec(Ep_up,Ep);
Ep_up = feval(PRF.M.IS, Ep_up, PRF.M, PRF.U, 'get_parameters');
Ep_up = spm_vec(Ep_up);
Ep_up = Ep_up(selected_idx);

% Latent variable lower CI -> underlying parameter
Ep_lo  = spm_vec(Ep);
Ep_lo(selected_idx) = ci_lower;
Ep_lo  = spm_unvec(Ep_lo,Ep);
Ep_lo = feval(PRF.M.IS, Ep_lo, PRF.M, PRF.U, 'get_parameters');
Ep_lo = spm_vec(Ep_lo);
Ep_lo = Ep_lo(selected_idx);

raw_string = sprintf('Raw param: %3.3f (95%% CI %3.3f to %3.3f)',...
    Ep_exp,Ep_lo,Ep_up);
set(handles.ci_text,'String',{ci_str;raw_string});

% Get predicted BOLD
y = feval(PRF.M.IS, Ep, PRF.M, PRF.U);

% Get residual
Y = PRF.Y;
R = Y.y(:,1) - y;
R = R - Y.X0*spm_inv(Y.X0'*Y.X0)*(Y.X0'*R);

% Plot BOLD
axes(handles.bold_axes);
cla;
x = [1:size(Y.y,1)]*Y.dt;
plot(x,y+R,'k');
hold on;
plot(x,y,'r');
title('BOLD Response','FontSize',16);
xlabel('Time (secs)'); ylabel('BOLD');

% List X,Y stimulus coordinates
pmax = PRF.M.pmax;
x_bins = -(pmax/2):0.1:(pmax/2);
y_bins = x_bins;
b      = length(x_bins);
[x2,y2] = meshgrid(x_bins,y_bins);
xy = [x2(:) y2(:)];

% Get PRF response
z_post   = feval(PRF.M.IS, Ep, PRF.M, PRF.U, 'get_response', xy);
z_post   = reshape(z_post,b,b);

% Plot PRF response
axes(handles.prf_axes);
cla;
h = imagesc(z_post);
axis square;
colormap jet;
colorbar;
if ~isempty(mask), set(h,'AlphaData',mask); end
title('PRF','FontSize',16);
set(gca,'YDir','normal','XTick',1:10:b,'YTick',1:10:b);
set(gca,'XTickLabel',x_bins(1:10:end),'YTickLabel',x_bins(1:10:end)); 
set(gca,'FontSize',12); 

axes(handles.ppd_axes);
cla;

if options.disp_ppd
    nsamp = options.nsamp;
    
    % Calculate predictive density
    pd = get_ppd(PRF,Ep,Cp,nsamp,xy);
    pd = reshape(pd,b,b);
    
    h = imagesc(pd);
    title('PRF with uncertainty','FontSize',16);
    axis square;
    colormap jet;
    colorbar;
    if ~isempty(mask), set(h,'AlphaData',mask); end
    set(gca,'YDir','normal','XTick',1:10:b,'YTick',1:10:b);
    set(gca,'XTickLabel',x_bins(1:10:end),'YTickLabel',x_bins(1:10:end)); 
    set(gca,'FontSize',12)
end

% -------------------------------------------------------------------------
function pd = get_ppd(PRF,Ep,Cp,nsamp,xy)

Ep_struct = Ep;

Ep = spm_vec(Ep);

pd = 0;

% Sample parameters
sEp = spm_normrnd(Ep, full(Cp), nsamp);

% f = gcf;
% figure;hist(sEp(3,:));
% figure(f);

% Run through likelihood
for i = 1:nsamp    
    % Sample from priors
    sEp_struct = spm_unvec(sEp(:,i),Ep_struct);
    
    g = feval(PRF.M.IS, sEp_struct, PRF.M, PRF.U, 'get_response', xy);
    pd = pd + g;
end

pd = pd .* (1/nsamp);

% -------------------------------------------------------------------------
function pselect_callback(varargin)
% Updates the display following a change of selected parameter

% TODO - get figure from control
f = gcf;

data = get(f,'UserData');

% Get newly selected parameter
handles   = data.handles;
new_param = get(handles.pselect,'Value');

% Store any changes
[Ep,Cp,options] = update_from_gui(f);
data.Ep = Ep;
data.Cp = Cp;
data.options = options;
data.sel_param_idx = new_param;
set(f,'UserData',data);

% Update controls
update_controls(data.Ep,data.Cp,data.options,data.handles);

% Display PRF
display_PRF();

% -------------------------------------------------------------------------
function slider_callback(varargin)

% TODO - get figure from control
f = gcf;

data = get(f,'UserData');
h    = data.handles;

% Transfer slider values to editors
set(h.exp_edit,'String',get(h.exp_slider,'Value'));
set(h.var_edit,'String',get(h.var_slider,'Value'));

display_PRF();
% -------------------------------------------------------------------------
function editor_callback(varargin)

% TODO - get figure from control
f = gcf;

data = get(f,'UserData');
h    = data.handles;

% Transfer editor values to sliders
set(h.exp_slider,'Value',str2double(get(h.exp_edit,'String')));
set(h.var_slider,'Value',str2double(get(h.var_edit,'String')));

display_PRF();