function spm_prf_review_surface(PRF,pidx,img)
% Displays PRF results on a cortical surface
%
% PRF  - estimated PRF structure
% pidx - index of the parameter to display within PRF.Ep{1}
% img  - (optional) displays this parameter map (.nii)

% Default hemisphere
hemi = 'lh';

alpha = 0.8;

if nargin < 3
    img = '';
end

if ~isempty(img)
    % There's no parameter index if an image was provided
    pidx = [];
end

% Select PRF dir
% -------------------------------------------------------------------------
if isfield(PRF,'dir') && ~isempty(PRF.dir)
    prf_dir = PRF.dir;
else
    [prf_dir,sts] = ...
        spm_select(1,'dir','Please select directory with lh_Srf.mat or rh_Srf.mat');
    if ~sts, return; end
end

% Assume PRF is in the GLM dir
glm_dir = prf_dir;

% Select vertices-voxels mapping files (created by spm_prf_import_surface)
% -------------------------------------------------------------------------
lh_coords = fullfile(prf_dir,'lh_Srf.mat');
rh_coords = fullfile(prf_dir,'rh_Srf.mat');

has_lh = exist(lh_coords,'file');
has_rh = exist(rh_coords,'file');

if ~has_lh && ~has_rh
    error('Could not find lh_Srf.mat or rh_Srf.mat in PRF.dir');
end

% Load space defining image
V = spm_vol(fullfile(glm_dir,'mask.nii'));

% Select parameter to display
% -------------------------------------------------------------------------
if isempty(pidx)
    param = '';
else
    pnames = fieldnames(PRF.Ep{1});
    param  = pnames{pidx};
end

% Create figure
% -------------------------------------------------------------------------
if isempty(img)
    title_str = sprintf('Surface: Parameter %s',param);
else
    title_str = 'Surface';
end
figure('Color','w','Name',title_str,'NumberTitle','off');

% Prepare text for drop-downs
% -------------------------------------------------------------------------
hemi_txt = {'Change hemisphere'};
if has_lh
    hemi_txt{end+1} = 'Left';
end
if has_rh
    hemi_txt{end+1} = 'Right';
end

export_txt = {'Export...','SamSrf'};

% Go!
% -------------------------------------------------------------------------
prepare_figure();
Cmap = display_brain();
display_key(Cmap);


function prepare_figure()
    % Adds controls to the figure
    spm_clf;

    uicontrol('Style','Popupmenu','Units','normalized', ...
        'Position',[0.01 0.88 0.30 0.1], 'String',hemi_txt,'FontSize',12,...
        'Callback',@hemisphere_callback);

    uicontrol('Style','Popupmenu','Units','normalized',  ...
        'Position',[0.01 0.80 0.30 0.1],'String',export_txt,'FontSize',12,...
        'Callback',@export_callback);
end

function Cmap = display_brain    
    % Load surface (created by spm_prf_import_surface.m)
    Srf = load(fullfile(glm_dir, [hemi '_Srf.mat']));
    Srf = Srf.Srf;
    
    % Unpack
    Vertices = Srf.Sphere;
    Faces    = Srf.Faces;
    Curv     = Srf.Curvature;
    fXYZ     = Srf.Voxels;
    
    % Clean up vertices
    Vertices(isnan(Vertices(:,1)),:) = NaN; 
    nvertex = length(Vertices);

    % Get voxel coordinates for pRF map
    pXYZmm = PRF.xY.XYZmm;
    nvox   = size(pXYZmm,2);
    pXYZ   = V.mat \ [pXYZmm;ones(1,length(pXYZmm))];
    pXYZ   = round(pXYZ(1:3,:))';

    % Prepare parameter map
    % ---------------------------------------------------------------------

    Y = zeros(1,nvox);
    Pp = zeros(1,nvox);    
    
    if isempty(img)
        % Get PRF parameters
        for i = 1:nvox
            p = feval(PRF.M.IS, PRF.Ep{i}, PRF.M, PRF.U, 'get_parameters');
            Y(i) = p.(param);
            [~,Pp(i)] = feval(PRF.M.IS, PRF.Ep, PRF.M, PRF.U, 'is_above_threshold', PRF.Cp, i, alpha);
        end
    else
        % Get parameters from provided nifti image
        
        % PRF mm -> vox
        prf_voxel_ind = sub2ind(V.dim,pXYZ(:,1),pXYZ(:,2),pXYZ(:,3));
        
        % Read parameter map
        Y = spm_read_vols(spm_vol(img));
        Y = Y(prf_voxel_ind);
        Pp = (Y ~= 0);    
        
        % Temporarily store in PRF, to support export to Srf callback
        PRF.img_Y = Y;
    end

    % Set pRF parameter for each vertex
    X  = zeros(nvertex,1);
    for i = 1:nvox

        % Identify vertices for this voxel
        vox_xyz = pXYZ(i,:);    
        tf  = (fXYZ(:,1) == vox_xyz(1)) & (fXYZ(:,2) == vox_xyz(2)) & ...
              (fXYZ(:,3) == vox_xyz(3));

        % Threshold parameter
        if Pp(i) > alpha
            X(tf) = Y(i);
        end
    end

    % Prepare surface curvature map
    % ---------------------------------------------------------------------

    % Greyscale for curvature
    CurvGrey = gray(11);         % Grey scale colour map
    CurvGrey = CurvGrey(1:10,:); % Remove black & white

    % Transform curvature
    Curv = -Curv + 0.5;
    Curv(Curv <= 0) = 0.000001;
    Curv(Curv > 1) = 1;
    Curv = ceil(Curv * size(CurvGrey,1));
    if length(Curv) == 1
        Curv = repmat(Curv, 1, size(Vertices,1));
    end

    % Set which parameters to include (remove very small ones)
    is_zero    = (abs(X) < 0.001);

    % Set colour of each vertex
    % ---------------------------------------------------------------------
    % Set bounds
    switch param
        case 'angle'
            current_min = -pi;
            current_max = pi;
            Cmap     = [circshift(fscol(100),75); CurvGrey];
        case 'dist'
            current_min = PRF.M.pmin;
            current_max = PRF.M.pmax/2;
            Cmap     = [flipud(jet(100)); CurvGrey];
        otherwise
            current_min = min(X(~is_zero));
            current_max = max(X(~is_zero));
            Cmap     = [jet(100); CurvGrey];
    end

    % Scale parameters to indices in the colormap (1 to 100)
    new_max = 100; new_min = 1;
    scale_factor = (current_max - current_min) / (new_max - new_min);
    X_scaled     = round(new_min + (X(~is_zero) - current_min) / scale_factor); 

    % Colours for each vertex (indices 100-110 are grayscale)
    Colours  = nan(nvertex,1);
    Colours(~is_zero) = X_scaled;
    Colours(is_zero)  = 100 + Curv(is_zero);

    % Display the mesh
    % ---------------------------------------------------------------------
    patch('vertices', Vertices, 'faces', Faces(:,[1 3 2]), ...
        'FaceVertexCData', Colours, 'FaceColor', 'interp', 'EdgeColor', 'none');
    axis off;

    % Set the view
    if strcmp(hemi,'lh')
        CamView = [42 -2 0.8];
    else
        CamView = [-50 -6 0.6];
    end
    set(gca, 'view', CamView(1:2));
    zoom(CamView(3));
    set(gca, 'projection', 'perspective');
    daspect([1 1 1]); % Correct aspect ratio

    colormap(Cmap);
    caxis([1 110]);

end

function display_key(Cmap)
    % Add colour wheel / colour bar
    xs = -10.5:0.05:10.5;
    ys = xs;

    np = length(xs);

    % x,y on-screen coordinates
    [x2,y2] = meshgrid(xs,ys);
    xy = [x2(:) y2(:)];

    % polar coordinates for each pixel
    dist  = sqrt( xy(:,1).^2 + xy(:,2).^2 );
    angle = atan2(xy(:,2),xy(:,1));

    % create images
    dist_map  = reshape(dist,np,np);
    angle_map = reshape(-angle,np,np);

    % threshold to stimulated area
    dist_map(dist_map < -10.5 | dist_map > 10.5) = nan;
    mask = ~isnan(dist_map);

    % background image
    bg = zeros(np,np);

    % create axes and plot
    ax=axes();
    set(ax,'Position',[0.01 0.3 0.2 0.2]);
    switch param
        case 'angle'
            h = imagesc(angle_map); 
            set(h,'AlphaData',mask);
        case 'dist'
            h = imagesc(dist_map); 
            set(h,'AlphaData',mask);
        otherwise
            colorbar;
    end
    axis off;
    axis square;
    set(gca,'XTick',[],'YTick',[]);
    colormap(ax,Cmap(1:100,:));
end

function hemisphere_callback(varargin)
    % Callback for the change hemisphere dropdown
    selected_str = varargin{1}.String(varargin{1}.Value);
    switch selected_str{1}
        case 'Left'
            hemi = 'lh';
        case 'Right'
            hemi = 'rh';
    end

    prepare_figure();
    Cmap = display_brain();
    display_key(Cmap);

end

function export_callback(varargin)
    % Callback from the export dropdown
    option = varargin{1}.String(varargin{1}.Value);
    
    if strcmp(option,'SamSrf')
        export_samsrf();
    end
end

function export_samsrf(varargin)
    % Exports the PRF results to Samsrf format    
    name = inputdlg('Please select a name','name',1);
    spm_prf_export_samsrf(PRF, hemi, name{1});
    disp('Export complete');
end

end