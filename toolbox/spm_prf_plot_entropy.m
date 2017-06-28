function fname = spm_prf_plot_entropy(PRF,fields,name,doplot)
% Computes and optionally displays a negative entropy map. More positive
% values (in units of nats) correspond to more confident parameter
% estimates.
%
% PRF    - estimated PRF structure
% fields - cell array of field names to include (default: 'all')
% name   - the filename will be 'entropy_name.nii' (default: 'entropy')
% doplot - whether to display the result (default: false)
%
% fname  - filename of generated entropy map
%
% ---------------------------------------------------------------------
% Copyright (C) 2017 Peter Zeidman
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

% Threshold for initial feature selection
alpha = 0.95;

if nargin < 2
    fields = 'all';
elseif ~iscell(fields)
    fields = cellstr(fields);
end

if nargin < 3
   name = 'entropy';
else
   name = ['entropy_' name];
end

if nargin < 4
    doplot = false;
end

nvox = length(PRF.Ep);

% Calculate
entropy = zeros(1,nvox); 
for v = 1:nvox    
   Cp = PRF.Cp{v};  
   Ep = PRF.Ep{v};
      
   idx = spm_fieldindices(Ep,fields{:});
   
   if isempty(idx)
       error('Field(s) not found');
   end
   
   Cp = Cp(idx,idx);
   
   % Feature selection
   [~,Pp] = feval(PRF.M.IS, PRF.Ep, PRF.M, PRF.U, 'is_above_threshold', PRF.Cp, v, alpha);
   
   if Pp > alpha
       entropy(v) = -spm_logdet(Cp);      
   end
end

% Create image
try 
    out_dir = PRF.dir;
catch
    out_dir = pwd;
end

SPM   = load_spm(PRF);
fname = create_results_img(PRF,entropy,name,SPM,out_dir);

% Display
if doplot
    spm_prf_review(PRF,'view_comparison',fname);
end

%--------------------------------------------------------------------------
function fname = create_results_img(PRF,entropy,image_name,SPM,dir_out)

out_file = [image_name '.nii'];
fname    = fullfile(dir_out,out_file);

if ~isfield(PRF,'Ep')
    error('Please ensure all models are estimated before continuing');
end

disp('Creating results image');

V = SPM.Vbeta(1);

% mm -> vox
Q   = ones(1,size(PRF.xY.XYZmm,2));
XYZ = round(V.mat \ [PRF.xY.XYZmm; Q]);

% vox -> idx
j = sub2ind(V.dim,XYZ(1,:),XYZ(2,:),XYZ(3,:));
j = round(j);

start_dir = pwd;
cd(dir_out);

% Make images
Y       = nan(V.dim);
Y(j)    = entropy;
V.n(1)  = 1;
V.fname = out_file;
spm_write_vol(V,Y);

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
catch
    [P,sts] = spm_select(1,'mat','Please select SPM.mat',{},pwd,'SPM.mat');
    if sts
        SPM = load(P);
    else
        return;
    end
end    
SPM = SPM.SPM;