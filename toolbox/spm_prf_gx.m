function [y] = spm_prf_gx(x,u,P,M)
% Simulated BOLD response to input
% FORMAT [y] = gx_hdm (x,u,P,M)
% y          - BOLD response (%)
% x          - state vector     (see my_fx_hdm)
% P          - Parameter vector 
% M          - model specification structure (see spm_nlsi)
%__________________________________________________________________________
%
% This function implements the BOLD signal model described in: 
%
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%
% Field-strength-specific parameters from:
%
% Heinzle J, Koopmans PJ, den Ouden HE, Raman S, and Stephan, KE (2016) 
% A hemodynamic model for layered BOLD signals. Neuroimage 125: 556-570.
%__________________________________________________________________________
%
% Adapted from spm_gx_hdm
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Karl Friston
 
% get MR field strength
%--------------------------------------------------------------------------
try
    B0 = M.B0;
catch
    B0 = 1.5;
end

% Biophysical constants
%==========================================================================
 
% time to echo (TE) (default 0.04 sec)
%--------------------------------------------------------------------------
try, TE = M.TE; catch, TE = 0.04; end
 
% resting venous volume (%)
%--------------------------------------------------------------------------
V0  = 4;

%P=spm_unvec(P,M.pE);

% estimated region-specific ratios of intra- to extra-vascular signal 
%--------------------------------------------------------------------------
ep  = 1*exp(P.epsilon);

% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation S:  R_iv = r0*[(1 - S)-(1 - S0)] (Hz)
%--------------------------------------------------------------------------
switch B0
    case 3
        r0 = 110;
    case 7
        r0 = 340;
    otherwise
        r0  = 25; % Backward compatibility
end
 
% frequency offset at the outer surface of magnetized vessels (Hz)
%--------------------------------------------------------------------------
switch B0
    case 1.5
        nu0 = 40.3; % Backward compatibility
    otherwise
        nu0 = 28.265 * B0;
end
 
% resting oxygen extraction fraction
%--------------------------------------------------------------------------
E0  = 0.4;
 
%-Coefficients in BOLD signal model
%==========================================================================
k1  = 4.3*nu0*E0*TE;
k2  = ep*r0*E0*TE;
k3  = 1 - ep;
 
%-Output equation of BOLD signal model
%==========================================================================
v   = exp(x(3));
q   = exp(x(4));
y   = V0*(k1.*(1 - q) + k2.*(1 - q./v) + k3.*(1 - v));
