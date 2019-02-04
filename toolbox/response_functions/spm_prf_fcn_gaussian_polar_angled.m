function varargout = spm_prf_fcn_gaussian_polar_angled(P,M,U,varargin)
% Gaussian PRF model with polar coordinates with hard constraints
%
% D - Degrees of visual angle used as units
% C - Constrained PRF width
%
% -------------------------------------------------------------------------
% FORMAT [y,Z] = spm_prf_fcn_gaussian(P,M,U)
% Return the BOLD and neuronal predicted timeseries
%
% P         parameters
% M,U       model, inputs
%
% y         fMRI time series
% Z         neuronal response
% -------------------------------------------------------------------------
% FORMAT P = spm_prf_fcn_gaussian(P,M,U,'get_parameters')
% Return the given parameters corrected for display
%
% P         parameters
% M,U       model, inputs
% -------------------------------------------------------------------------
% FORMAT S = spm_prf_fcn_template(P,M,U,'get_summary')
% Summarises the pRF with simple (Gaussian) parameters x,y,width,beta
%
% S         structure with fields x,y,width,beta
% M,U       model, inputs
% -------------------------------------------------------------------------
% FORMAT P = spm_prf_fcn_gaussian(P,M,U,'get_response',xy)
% Return the instantaneous response of the PRF at coordinates xy
%
% P         parameters
% M,U       model, inputs
% xy        [2xn] vector of coordinates to evaluate the PRF
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

if nargin < 4
    P = correct_parameters(P,M);
    [z,varargout{1},varargout{2}] = integrate_model(P,M,U);
else
    action = varargin{1};

    switch action         
        case 'get_parameters'
            % Gets the parameters for display
            
            % Correct neuronal parameters
            P = correct_parameters(P,M);
            
            % Corrected haemodynamic parameters
            P.transit = exp(P.transit);
            P.decay   = exp(P.decay);
            P.epsilon = exp(P.epsilon);                  
            
            varargout{1} = P;
        case 'get_summary'
            % Gets x,y,sigma,beta
            
            P = correct_parameters(P,M);
            
            x     = P.dist .* cos(P.angle);
            y     = P.dist .* sin(P.angle);            
            width = max(P.width_x, P.width_y);
            beta  = P.beta;
            
            varargout{1} = struct('x',x,'y',y,'width',width,'beta',beta);            
        case 'is_above_threshold'
            Cp = varargin{2};
            v  = varargin{3};
            if length(varargin) >= 4
                alpha = varargin{4};
            else
                alpha = 0.95;
            end
            varargout{1} = is_above_threshold(P,M,Cp,v,alpha);
        case 'get_response'
            xy = varargin{2};                        
            P = correct_parameters(P,M);
            varargout{1} = get_response(P,M,U,xy,false);
        case 'get_priors'
            [pE,pC] = get_priors(M);
            varargout{1} = pE;
            varargout{2} = pC;
        case 'glm_initialize'
            y = varargin{2};
            pE = glm_initialize(P,M,U,y);
            varargout{1} = pE;
        otherwise
            error('Unknown action');
    end
end

% -------------------------------------------------------------------------
function [z,y,Z] = integrate_model(P,M,U)
% Integrates the neuronal and BOLD models

n  = length(U);            % Number of volumes = inputs
%T0 = M.T0;                 % Microtime offset

z = zeros(1,U(1).nbins);

for t = 1:n        
    % Microtime index for this volume
    %ind = U(t).ind(T0);
    ind = U(t).ind;

    % Get coordinates of all stimulated pixels
    xy = [U(t).dist U(t).angle];
   
    % Get PRF response
    x = get_response(P,M,U,xy,true);
   
    % Sum over pixels
    z(ind) = z(ind) + sum(x);   
end

if nargout > 1
    % Integrate BOLD model
    Z.u=z';
    Z.dt=M.dt;
    y=spm_int(P,M,Z);
end

% -------------------------------------------------------------------------
function x = get_response(P,M,U,xy,is_polar)
% Get PRF response to a stimulus covering the given pixels
%
% P  - parameters
% xy - [2xn] vector of pixel coordinates
% is_polar - true if the xy coordinates are polar

if isempty(xy)
    x = 0;
    return;
end

% Inputs: polar coords -> x,y
if is_polar
    dist    = xy(:,1);
    angle   = xy(:,2);
    xy(:,1) = dist .* cos(angle); % x
    xy(:,2) = dist .* sin(angle); % y
end

% Parameters: polar coords -> x,y
mu_x = P.dist .* cos(P.angle);
mu_y = P.dist .* sin(P.angle);

% Std dev -> covariance matrix
sigma = diag( [P.width_x P.width_y] );
sigma = sigma .* sigma;

% Correlation (angle) centre
sigma_corr = P.rotation * P.width_x * P.width_y;
sigma      = sigma + [0 sigma_corr; ...
                      sigma_corr 0];

% Normalized Gaussian
x = spm_mvNpdf(xy', [mu_x mu_y], sigma);

% Scale
x = P.beta .* x;

% -------------------------------------------------------------------------
function x2 = constrain_parameter(x,d_min,d_max)
% Convert latent variable x to parameter x2 where d_min <= x <= d_max
% Gives uniform priors with mean 0 and variance 1
%
% x     - real-valued parameters
% d_min - target minimum value
% d_max - target maximum value
%
% x2 - parameters scaled between d_min and d_max

d_range = d_max - d_min;
x2      = (d_range .* spm_Ncdf(x,0,1)) + d_min;

% -------------------------------------------------------------------------
function P = correct_parameters(P,M)
% Prepare neuronal parameters for use e.g. exponentiate log parameters
%
% P - parameters
% M - model

% Constrain PRF centre
radius    = M.pmax / 2;
P.dist    = constrain_parameter(P.dist, M.rmin, radius);
P.angle   = constrain_parameter(P.angle, -pi, pi);

% Convert log beta -> beta
P.beta    = exp(P.beta);

% Constrain rotation
P.rotation  = constrain_parameter(P.rotation,-1,1);

% Scale width (SD) between 0.1 degrees and half stimulus diameter
P.width_x = constrain_parameter(P.width_x, M.pmin, M.pmax / 2);
P.width_y = constrain_parameter(P.width_y, M.pmin, M.pmax / 2);

% -------------------------------------------------------------------------
function tf = is_above_threshold(P,M,Cp,v,alpha)
% Evaluate whether the model with parameters P and covariance Cp passes an
% arbitrary threshold for display
%
% P     - parameters
% M     - model
% Cp    - covariance matrix
% alpha - threshold

pC = M.pC{v};
pE = M.pE{v};
Cp = Cp{v};
Ep = P{v};

rE = spm_vec(pE);
rC = spm_vec(pC);

% Indices of parameters to keep in nested model
np = length(rE);
q  = zeros(1,np);
q(end-2:end) = 1; % haemo
q(5)        =  1; % beta

% Remove others
rC(q~=1)=0;

% BMR
F = spm_log_evidence(Ep,Cp,pE,diag(spm_vec(pC)),rE,diag(rC));

% Convert to post. prob - from spm_dcm_compare
F    = [F 0];
i    = F < (max(F) - 32);
P    = F;
P(i) = max(F) - 32;
P    = P - min(P);
P    = exp(P);
P    = P/sum(P);

tf = P(2) > alpha;

% -------------------------------------------------------------------------
function [pE,pC] = get_priors(M)

% Get the neuronal priors
pE.dist    = 0;         pC.dist    = 1;
pE.angle   = 0;         pC.angle   = 1;
pE.width_x = 0;         pC.width_x = 1;
pE.width_y = 0;         pC.width_y = 1;
pE.beta    = -2;        pC.beta    = 5;
pE.rotation= 0;         pC.rotation= 1;

% -------------------------------------------------------------------------
function P = glm_initialize(P,M,U,y)

% Coordinates as latent variable dist / angle
x     = -2:0.2:2;
sigma = P.width_x;

% Candidates for mu
[dist, angle] = meshgrid(x,x);
mu = [dist(:) angle(:)];

% Model space
n_sigma = size(sigma,2);
n_mu    = size(mu,1);
sigma   = kron(ones(n_mu,1), sigma);
sigma   = sigma(:);
k       = [repmat(mu,n_sigma,1) sigma];

nm = size(k,1);
sse = zeros(1, nm);

for i = 1:nm
    
    P2         = P;
    P2.dist    = k(i,1);
    P2.angle   = k(i,2);
    P2.width_x = k(i,3);
    P2.width_y = k(i,3);
    P2.beta    = 0;
    
    P2 = correct_parameters(P2,M);
    
    %fprintf('M%d  dist: %2.2f angle %2.2f\n',i,P2.dist,P2.angle);
    
    % Get neuronal response
    z = integrate_model(P2,M,U);
    
    % Convolve with canonical HRF
    XU.u=z';
    XU.dur=0;
    XU.dt=U(1).dt;
    XU.name{1} = 'gPRF';
    xBF.name   = 'hrf';
    xBF.order  = 1;
    xBF.length = 32;
    xBF.dt     = U(1).dt;
    xBF        = spm_get_bf(xBF);    
    z = spm_Volterra(XU, xBF.bf, 1);
    
    % Downsample regressor
    fMRI_T = spm_get_defaults('stats.fmri.t');
    ns     = M.ns;
    ind    = (0:(ns - 1))*fMRI_T + M.T0;
    z      = z(ind,:);
    
    X = [z ones(ns,1)];
    
    % Invert
    beta   = pinv(X) * y;
    e      = y - X*beta;
    sse(i) = e'*e;   
end

[i,idx] = min(sse);

P.dist    = k(idx,1);
P.angle    = k(idx,2);
%P.width   = k(idx,3);