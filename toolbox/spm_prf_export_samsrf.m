function Srf = spm_prf_export_samsrf(PRF, hemi, name)
% Exports an estimated PRF to SamSrf format
%
% PRF  - PRF structure (estimated)
% hemi - hemisphere ('lh' or 'rh')
% name - name for the file
%
% Outputs:
%
% Srf.Values  - field names of estimated parameters {p x 1}
% Srf.Data    - vertex-wise estimated parameters [p x v]
% Srf.Roi     - vertex index for each vertex in Srf.Data [v x 1]
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

if nargin < 3 || isempty(name)
    name = 'BayespRF';
end

% Load surface data (created by spm_prf_import_surface.m)
Srf = load(fullfile(PRF.dir,[hemi '_srf.mat']));
Srf = Srf.Srf;

% Parameters
Srf.Values = {'R^2','x0','y0','Sigma','Beta','Field sign'}';

nvox  = length(PRF.Ep);
nvert = length(Srf.Vertices);

% Get voxel coordinates for each voxel in the PRF
XYZmm = PRF.xY.XYZmm;
M     = PRF.xY.spec.mat;
XYZ   = zeros(4,size(XYZmm,2));
for v = 1:nvox
    XYZ(:,v) = round(M \ [XYZmm(:,v);1]);
end
XYZ = XYZ(1:3,:);

Data = []; % Estimated parameters
Roi  = []; % Indices of relevant vertices

disp('Exporting');
disp('Note: The posterior probability of the model will be used in place of R^2');

% Prepare progress bar and bring to front
spm_progress_bar('Init',nvox,'','Voxel');
Finter = spm_figure('FindWin','Interactive');
if ~isempty(Finter)
    figure(Finter);
end

for v = 1:nvox
    
    spm_progress_bar('Set',v);
    
    % Identify vertices corresponding to this voxel
    % ---------------------------------------------------------------------
    
    % Compute distance between each PRF XYZ and each vertex XYZ
    d = Srf.Voxels - repmat(XYZ(:,v)',nvert,1);
    
    % Identify vertices which match
    j = (d(:,1) == 0 & d(:,2) == 0 & d(:,3) == 0);
    
    vertex_indices = find(j);
    
    % Store
    Roi = [Roi; vertex_indices];    

    % Get parameters
    % ---------------------------------------------------------------------
    
    % Get x, y, width, beta
    S = feval(PRF.M.IS, PRF.Ep{v}, PRF.M, PRF.U, 'get_summary');
    
    % Calculate R^2 (slow!)
    %predicted = feval(PRF.M.IS, PRF.Ep{v}, PRF.M, PRF.U);
    %r2 = corrcoef(predicted, PRF.Y.y(:,v)) .^ 2;
    %r2 = r2(1,2);
    
    % In place of the R^2, use the model probability
    [~, Pp] = feval(PRF.M.IS, PRF.Ep, PRF.M, PRF.U, 'is_above_threshold', PRF.Cp, v);
    r2 = Pp;    
    
    % What's this then?
    field_sign = 1;
    
    % Assemble parameter vector
    r = [r2 S.x S.y S.width S.beta field_sign]';
    
    % Append optional data in PRF.img_Y (added by spm_prf_review_surface.m
    % when viewing a parameter map, e.g. model comparison results)
    if isfield(PRF,'img_Y')
        r = [r; PRF.img_Y(v)];
    end
    
    % Replicate vector for each relevant vertex
    Data = [Data repmat(r,1,length(vertex_indices))];
    
end

if isfield(PRF,'img_Y')
    Srf.Values{end+1} = name;
end

spm_progress_bar('Clear');

% Pack
Srf.Data = Data;
Srf.Hemisphere = hemi;
Srf.Roi = Roi;

% Maximum stimulus radius
Eccentricity = (PRF.M.pmax / 2); %#ok<NASGU>

% Save
fn = sprintf('%s_Srf_%s.mat',hemi,name);
fn = fullfile(PRF.dir,fn);
save(fn,'Srf','Eccentricity');
fprintf('Saved %s\n',fn);