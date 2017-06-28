function varargout = spm_prf_fcn_template(P,M,U,varargin)
% Template pRF response function. This example models neuronal
% activity as a scaled version of the inputs:
%
% z(t) = alpha * u(t)
%
% Where alpha is a parameter estimated from the data.
% 
%
% Inputs:
%
% P      - parameter structure
% M      - model structure
% U      - experimental timing
% action - (optional) the action to performm. If not provided, the
%          predicted BOLD and neuronal timeseries are returned.
%
% -------------------------------------------------------------------------
% FORMAT [y,Z] = spm_prf_fcn_template(P,M,U)
% Return the BOLD and neuronal predicted timeseries
%
% P         parameters
% M,U       model, inputs
%
% y         fMRI time series
% Z         neuronal response
% -------------------------------------------------------------------------
% FORMAT P = spm_prf_fcn_template(P,M,U,'get_parameters')
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
% FORMAT tf = spm_prf_fcn_template(P,M,U,'is_above_threshold',Cp,v)
% Return whether the model with parameters P and covariance Cp passes an
% arbitrary threshold for display
%
% P         parameters
% M,U       model, inputs
% Cp        parameter covariance matrix
% v         voxel index
% -------------------------------------------------------------------------
% FORMAT x = spm_prf_fcn_template(P,M,U,'get_response',xy)
% Return the instantaneous response of the PRF at coordinates xy
%
% P         parameters
% M,U       model, inputs
% xy        [2xn] vector of coordinates to evaluate the PRF
% -------------------------------------------------------------------------
% FORMAT [pE,pC] = spm_prf_fcn_template(P,M,U,'get_priors')
% Return the priors for the model. Importantly, this defines which
% parameters are in the model.
%
% pE        structure or vector of prior expectations
% pC        prior covariance maitrx
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
    % Integrate the model over time and return neuronal timeseries z and 
    % BOLD timeseries y. As an example, here we have neuronal model 
    % z(t) = alpha, where alpha is an estimated parameter.
    
    % Number of volumes = inputs
    n  = length(U); 

    % Neural timeseries. A vector with one entry per microtime bin. The 
    % U.nbins field is injected automatially by spm_prf_analyse
    z = zeros(1,U(1).nbins);
        
    for t = 1:n    
        % Microtime index for this volume
        ind = U(t).ind;        
        
        % pRF response
        z(ind) = P.alpha;
    end
            
    % Integrate the BOLD model
    Z.u=z';
    Z.dt=M.dt;
    y=spm_int(P,M,Z);      
    
    varargout{1} = y;
    varargout{2} = Z;
else
    % This section of the code provides information on the model, primarily
    % for plotting purposes in spm_prf_review()
    
    action = varargin{1};

    switch action         
        case 'get_parameters'
            % Get the parameters with any corrections needed for 
            % display            
            varargout{1} = P;
        case 'get_summary'
            % Get a summary of the pRF shape under Gaussian assumptions
            varargout{1} = ...
                struct('x',P.x,'y',P.y,'width',P.width,'beta',P.beta);
        case 'is_above_threshold'
            % Return binary vector identifying whether each voxel is
            % above some threshold for display            
            varargout{1} = 1;
        case 'get_response'            
            % Return the prediction of the model at coordinates xy            
            xy = varargin{2};            
            varargout{1} = ones(1,length(xy));            
        case 'get_priors'
            % Return a structure containing the priors
            pE.alpha = 1;
            pC.alpha = 1;
            varargout{1} = pE;
            varargout{2} = pC;            
        case 'glm_initialize'                        
            % (Optional) Return parameters initialized using some
            % rapid initial search on timeseries y
            y = varargin{2};
            
            varargout = P; 
        otherwise
            error('Unknown action');
    end
end