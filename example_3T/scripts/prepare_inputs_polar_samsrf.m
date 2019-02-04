function U = prepare_inputs_polar_samsrf(ApFrm,TR, nmicrotime, stim_duration, stim_diameter, varargin)
% Produces the input structure needed for spm_prf_analyse() given a 3D
% stimuli matrix with dimensions [x,y,t].
%
% Inputs:
%
% ApFrm         - [x,y,t] binary matrix indiciating which pixel locations
%                 were illuminated at each time point t
% TR            - Scanner repetition time (TR)
% nmicrotime    - Bins per TR
% stim_duration - Duration of stimuli (secs)
% stim_diameter - Diameter of stimuli in degrees
% pmin          - (optional) Minimum pRF size to be modelled, default 0.5
% rmin          - (optional) Minimum radius of stimulated area, default 0
% imscale       - (optional) Scale of image space for modelling, default 41

%
% Returns:
%
% U           - Input structure to feed to spm_prf_analyse.m

% Settings
nmicrotime    = 16;     % Bins per TR
stim_duration = 1;      % Duration of stimuli (secs)
stim_diameter = 17;     % Diameter of stimuli in degrees

n = size(ApFrm,3);

% Cechk for optional inputs
numvarargs = length(varargin);
if numvarargs > 1
    error('varargin:TooManyInputs', ...
          'requires at most 1 optional inputs');
end

% set defaults for optional inputs
optargs = {0.5, 0, 41};

% now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[pmin, rmin, imscale] = optargs{:};

for t = 1:n

    % Read
    im = ApFrm(:,:,t);

    % Skip if a baseline volume
    if ~any(im)
        continue;
    end

    % Rescale to imscale x imscale resolution
    res = [imscale imscale];
    im  = imresize(im,res);

    % Binarize
    im = im > 0.01;

    % Extract stimulated coordinates
    [y,x] = ind2sub(res,find(im));

    % Rescale [1,imscale] to units of visual angle (degrees)
    r = (stim_diameter/2);
    x = rescale(x, 1, res(1), -r, r);
    y = rescale(y, 1, res(1), -r, r);

    % Flip the y axis so positive is up
    y = -y;

    % Convert from x,y to distance >= 0 and angle (-pi,pi]
    dist  = sqrt( (x .^ 2) + (y .^ 2) );
    angle = atan2(y,x);

    U(t).dist  = dist;         % Distance
    U(t).angle = angle;        % Angle
    U(t).ons = TR * (t-1);     % Onset (secs)
    U(t).dur = stim_duration;  % Duration (secs)
    U(t).dt  = TR/nmicrotime;  % 1 bin=1/16 second, i.e. 16 bins per second
    U(t).pmax = stim_diameter; % Stimulus diameter
    U(t).pmin = pmin;          % Minimum PRF size
    U(t).rmin = rmin;          % Minimum radius of stimulated area

    t = t + 1;
end

% -------------------------------------------------------------------------
function new_X = rescale(X, current_min, current_max, new_min, new_max)
% Rescale X to the range [new_min new_max]
scale_factor = (current_max - current_min) / (new_max - new_min);
new_X = new_min + (X - current_min) / scale_factor;
