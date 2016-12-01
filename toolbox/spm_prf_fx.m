function [y] = spm_prf_fx(x,u,P,M)
% State equation for Neuro- and Hemo-dynamics
% FORMAT [y] = fx_hdm (x,u,P,M)

%   x - state vector
%   u(i) - drive to ith population
%
%   x(1) - vasodilatory signal                   s
%   x(2) - rCBF                                  ln(f)
%   x(3) - venous volume                         ln(v)
%   x(4) - deoyxHb                               ln(q)
%   x(5:end) - neuronal activity for each condition
%
%   y      - dx/dt
%
%___________________________________________________________________________
%
% References for hemodynamic & neuronal state equations:
% 1. Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
%    changes during brain activation: The Balloon model. MRM 39:855-864,
%    1998.
% 2. Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in
%    fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
%    Neuroimage 12:466-477, 2000.
%__________________________________________________________________________
%
% Adapted from spm_fx_hdm
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Karl Friston

% Hemodynamic motion
%==========================================================================

% hemodynamic parameters
%--------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) - autoregulation                                 d(ds/dt)/df)
%   H(3) - transit time                                   (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal
%--------------------------------------------------------------------------

H        = [0.64 0.32 2.00 0.32 0.32];

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(2:4) = exp(x(2:4));

% signal decay
%--------------------------------------------------------------------------
sd       = H(1)*exp(P.decay);

% transit time
%--------------------------------------------------------------------------
tt       = H(3)*exp(P.transit);

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv       = x(3).^(1/H(4));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff       = (1 - (1 - H(5)).^(1./x(2)))/H(5);

% Total neural activity
ztot = sum(u);

% Hemodynamics
%--------------------------------------------------------------------------
y(1)   = ztot - sd.*x(1) - H(2)*(x(2) - 1);
y(2)   = x(1)./x(2);
y(3)   = (x(2) - fv)./(tt.*x(3));
y(4)   = (ff.*x(2) - fv.*x(4)./x(3))./(tt.*x(4));
y      = y(:);