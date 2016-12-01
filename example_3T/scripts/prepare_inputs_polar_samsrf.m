function U = prepare_inputs_polar_samsrf(ApFrm,TR)
% Returns coordinates of pixels of the screen stimulated at each time point
%
% P           - [x,y,t] matrix
% TR          - Scanner repetition time (TR)
% num_session - Which session to use

% Settings
nmicrotime    = 16;     % Bins per TR
stim_duration = 1;      % Duration of stimuli (secs)
stim_diameter = 17;     % Diameter of stimuli in degrees

n = size(ApFrm,3);

for t = 1:n
    
    % Read
    im = ApFrm(:,:,t);
    
    % Skip if a baseline volume
    if ~any(im)
        continue;
    end
    
    % Rescale to 41 x 41 resolution
    res = [41 41];
    im  = imresize(im,res);
    
    % Binarize
    im = im > 0.01;
    
    % Extract stimulated coordinates
    [y,x] = ind2sub(res,find(im));
    
    % Rescale [1,41] to to units of visual angle (degrees)
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
    U(t).dt  = TR/nmicrotime;  % 1 bin=1/8 second, i.e. 8 bins per second
    U(t).pmax = stim_diameter; % Stimulus diameter
    U(t).pmin = 0.5;           % Minimum PRF size

    t = t + 1;
end

% -------------------------------------------------------------------------
function new_X = rescale(X, current_min, current_max, new_min, new_max)   
% Rescale X to the range [new_min new_max]
scale_factor = (current_max - current_min) / (new_max - new_min);
new_X = new_min + (X - current_min) / scale_factor;    