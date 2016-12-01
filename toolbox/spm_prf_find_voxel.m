function [idx,distmm] = spm_prf_find_voxel(XYZmm,PRF)
% Identifies the index of the voxel in the PRF with coordinates closest to
% those provided
%
% XYZmm - [3 x 1] vector of coordinates (mm)
% PRF   - PRF structure
%
% idx    - Index of the closest voxel
% distmm - Euclidean distance in mm of closest voxel to the target
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

XYZmm = XYZmm(:)';

mm = PRF.xY.XYZmm;

xy_dist = [];
for i = 1:size(mm,2)
    xy_dist(i) = pdist2(mm(:,i)',XYZmm);
end

[~,idx] = min(xy_dist);
distmm  = xy_dist(idx);