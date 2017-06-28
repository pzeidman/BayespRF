function [BPA,xY,XYZmm] = spm_prf_bpa_within_subject(PRF,xY,nocond)
% Bayesian Parameter Averaging over voxels.
%
% PRF    - subject's estimated pRF model
% xY     - Optional ROI definition for performing averaging 
%          (see spm_regions.m) or vector of voxel indices
% nocond - If true, disables conditional dependencies between parameters
%
% Returns:
% BPA    - pRF model averaged over voxels
% xY     - Updated ROI definition (optional)
% XYZmm  - Updated mm coordinates (optional)
%
% For details on why one might want to disable conditional dependencies,
% see https://en.wikibooks.org/wiki/SPM/Bayesian_Parameter_Averaging_(BPA)
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

if ischar(PRF)
    PRF = load(PRF);
    PRF = PRF.PRF;
end

nvox = length(PRF.Ep);

if nargin < 2
    % Include all voxels
    included_voxels = 1:nvox;
else
    % Mask
    if isnumeric(xY)
        included_voxels = xY;
    else
        if ischar(xY)
            xY = struct('def','mask','spec',xY);
        end
        [xY, XYZmm, included_voxels] = spm_ROI(xY, PRF.xY.XYZmm);
    end
end

if nargin < 3
    nocond = false;
end

% Extract posteriors and priors
nincluded = length(included_voxels);
GCM = cell(nincluded,1);
i = 1;
for v = included_voxels
    
    if v == included_voxels(1)
        GCM{i} = PRF;
    end
    
    GCM{i}.Cp = PRF.Cp{v};
    GCM{i}.Ep = PRF.Ep{v};
    
    GCM{i}.M.pC = PRF.M.pC{v};
    GCM{i}.M.pE = PRF.M.pE{v};
    i = i + 1;
end

% Average
BPA = spm_dcm_bpa(GCM,nocond);

% Use the cell array convention of PRF
BPA.Ep = {BPA.Ep};
BPA.Cp = {BPA.Cp};
BPA.M.pE = {BPA.M.pE};
BPA.M.pC = {BPA.M.pC};