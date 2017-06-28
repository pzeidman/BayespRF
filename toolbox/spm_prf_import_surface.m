function Srf = spm_prf_import_surface(glm_dir, struct_img, surf_dir, hemi)
% Imports Freesurfer surfaces to Matlab (SamSrf format)
%
% glm_dir        - directory containing SPM.mat
% struct_img     - filename of the subject's structural MRI
% surf_dir       - directory containing Freesurfer surfaces
% hemi           - 'lh' or 'rh' for left or right hemisphere respectively
%
% Outputs:
% Srf.Vertices   - from Freesurfer [v x 3]
% Srf.Curvature  - from Freesurfer [1 x v]
% Srf.Faces      - from Freesurfer [v x 3]
% Srf.Pial       - from Freesurfer 
% Srf.Inflated   - from Freesurfer 
% Srf.Sphere     - from Freesurfer 
% Srf.Voxels     - XYZmm of each vertex in EPI space [v x 3]
% Srf.Hemisphere - left or right hemisphere 'lh' or 'rh'
%
% The .mat file is saved in the GLM directory (<hemi>_Srf.mat), together 
% with a .nii image of voxels on the surface (<hemi>_surface.nii).
%
% Adapted from samsrf_vol2srf from the SamSrf toolbox by Sam Schwarzkopf
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

% Use location half way through the cortical layers
ctxstep = 0.5;
    
% Prepare image headers
% -------------------------------------------------------------------------

spm_mask_img = fullfile(glm_dir,'mask.nii');

% Read registration matrices from FreeSurfer
coreg = fullfile(surf_dir,'Coregistration.txt');
if exist(coreg, 'file')
    fs2nii = dlmread(coreg);
    Tmov   = fs2nii(1:4,:);
    Reg    = fs2nii(5:8,:);
    useRegDat = true;
else
    useRegDat = false;
end

% Load structural image
hdr = spm_vol(struct_img);

% Origin in the actual structural
nii_orig = hdr.mat(1:3,4);

% Origin in Freesurfer space (1/2 dimensions)
fs_orig = hdr.dim' / 2;
fs_orig = fs_orig([3 1 2]) .* sign(nii_orig); 

% Load functionals
fhdr = spm_vol(spm_mask_img);

% Adjust transformation matrix
% -------------------------------------------------------------------------
mov = nii_orig - fs_orig;
mat = fhdr.mat;
if useRegDat
    smat = hdr.mat;
else
    mat(1:3,4) = mat(1:3,4) - mov;
end

% Load surface vertices
% -------------------------------------------------------------------------
[V0,F] = fs_read_surf(fullfile(surf_dir, [hemi '.white']));      % Grey-white surface
P      = fs_read_surf(fullfile(surf_dir, [hemi '.pial']));       % Pial surface
I      = fs_read_surf(fullfile(surf_dir, [hemi '.inflated']));   % Inflated surface
S      = fs_read_surf(fullfile(surf_dir, [hemi '.sphere']));     % Spherical surface
C      = Read_FreeSurfer(fullfile(surf_dir,[hemi '.curv.asc'])); % Curvature 

N      = P - V0; % Normal vectors for each vertex 

% Get voxel coordinates
% -------------------------------------------------------------------------

% Step through cortex layers
V = V0 + N*ctxstep;

% Transform vertices (mm) to voxel coordinates
if useRegDat    
    % Inverse matrix of operations performed by Freesurfer
    M = Tmov \ Reg;
    
    % Vertex mm -> structural voxels
    sXYZ = M * [V'; ones(1,size(V,1))]; 
    
    % Structural voxels -> mm
    fXYZmm = smat * sXYZ;
    
    % mm -> functional voxels
    fXYZ = mat \ fXYZmm; 
           
    % Clean up
    fXYZ = round(fXYZ);
else
    % Functional mm -> vox
    fXYZ = round(mat \ [V'; ones(1,size(V,1))]); 
end
fXYZ = fXYZ(1:3,:)';

% Save in Samsrf format
% -------------------------------------------------------------------------
Srf = struct();
Srf.Voxels     = fXYZ;
Srf.Vertices   = V0;
Srf.Pial       = P;
Srf.Inflated   = I;
Srf.Sphere     = S;
Srf.Curvature  = C(:,5)';
Srf.Faces      = F;
Srf.Hemisphere = hemi;

% Stop here if the Srf is to be returned as a matrix
if nargout > 0 
    return
end

save(fullfile(glm_dir,[hemi '_Srf.mat']),'Srf');

% Create mask image
% -------------------------------------------------------------------------

% Remove out of range coordinates
is_ok = fXYZ(:,1)>0 & ...
        fXYZ(:,2)>0 & ...
        fXYZ(:,3)>0 & ...
        fXYZ(:,1)<fhdr(1).dim(1) & ...
        fXYZ(:,2)<fhdr(1).dim(2) & ...
        fXYZ(:,3)<fhdr(1).dim(3);        
fXYZ(~is_ok,:) = [];

Y = nan(fhdr.dim);
i = sub2ind(fhdr.dim,fXYZ(:,1),fXYZ(:,2),fXYZ(:,3));
Y(i) = 1;

start_dir = pwd;
cd(glm_dir);
fhdr.fname = [hemi '_surface.nii'];
fhdr = rmfield(fhdr,'pinfo');
spm_write_vol(fhdr,Y);
cd(start_dir);