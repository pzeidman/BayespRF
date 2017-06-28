function spm_prf_import_label(label_fn, glm_dir)
% Imports a Freesurfer label into Nifti format. 
% Requires an existing lh_Srf.mat or rh_srf.mat to exist in the subject's
% GLM directory (created by spm_prf_import_surface.m)
%
% label_fn - filename of the Freesurfer label
% glm_dir  - Directory containing surface .mat file
%
% Based on code from the SamSrf toolbox:
% http://dx.doi.org/10.6084/m9.figshare.1344765.v23
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

% Get labelname
[pathstr,label_name] = fileparts(label_fn);

% Choose hemisphere
hemi = label_name(1:2);
if ~strcmp(hemi,'lh') && ~strcmp(hemi,'rh')
    error('Label should begin with lh_ or rh_');
end

% Load Srf file (from spm_prf_import_surface.m)
srf_file = fullfile(glm_dir, sprintf('%s_Srf.mat',hemi));
Srf = load(srf_file);
Srf = Srf.Srf;

% Read vertex IDs from Freesurfer label
V = Read_FreeSurfer(label_fn);
V = V(:,1)+1;

% Identify coordinates matching vertices
XYZ = Srf.Voxels(V,:);

% Remove out of range voxels
spm_mask_img = fullfile(glm_dir,'mask.nii');
V = spm_vol(spm_mask_img);

is_ok = XYZ(:,1)>0 & ...
        XYZ(:,2)>0 & ...
        XYZ(:,3)>0 & ...
        XYZ(:,1)<V.dim(1) & ...
        XYZ(:,2)<V.dim(2) & ...
        XYZ(:,3)<V.dim(3);    
XYZ = XYZ(is_ok,:);

% Enter output dir
start_dir = pwd;
cd(glm_dir);

% Write image
V.fname = sprintf('%s.nii',label_name);

Y = nan(V.dim);
for i= 1:size(XYZ,1)
    Y( sub2ind(V.dim,XYZ(i,1),XYZ(i,2),XYZ(i,3)) ) = 1;
end
spm_write_vol(V,Y);

cd(start_dir);

end