function [spE,sEp,prior_pd,post_pd] = spm_prf_get_ppd(PRF,xy,idx,nsamp)
% Sample from the prior (or posterior) of a pRF model and run samples 
% through the likelihood, to give the prior (or posterior) predictive 
% density (PPD).
%
% PRF   - estimated PRF
% xy    - coordinates in stimulus space to sample
% idx   - index of the PRF within the PRF structure
% nsamp - the number of samples to draw
%
% Returns:
% spE      - samples from the prior [parameters x samples]
% sEp      - samples from the posterior [parameters x samples]
% prior_pd - average sampled PRF response (prior predictive density)
% prior_pd - average sampled PRF response (posterior predictive density)
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
    nsamp = 1;
end

% Sample from priors and posteriors
% -------------------------------------------------------------------------

pE_struct = PRF.M.pE{idx};

% Get priors and posteriors
pE = spm_vec(PRF.M.pE{idx});
pC = spm_vec(PRF.M.pC{idx});
Ep = spm_vec(PRF.Ep{idx});
Cp = full(PRF.Cp{idx});

if isvector(pC), pC = diag(pC); end

% Sample
spE = spm_normrnd(pE, pC, nsamp);
sEp = spm_normrnd(Ep, Cp, nsamp);

if nargout < 3
    return;
end

% Integrate model with sampled parameters
% -------------------------------------------------------------------------

% Create progress window
f = gcf;
f_progress = spm_figure('GetWin', 'Interactive');
spm_progress_bar('Init',nsamp,[],'Samples');

prior_pd = 0;
post_pd  = 0;

for i = 1:nsamp    
    
    if mod(nsamp,100) == 0
        spm_progress_bar('Set',i);
    end
    
    % Integrate model using prior sample
    spE_struct = spm_unvec(spE(:,i),pE_struct);
    
    g = feval(PRF.M.IS, spE_struct, PRF.M, PRF.U, 'get_response', xy);
    prior_pd = prior_pd + g;
    
    % Integrate model using posterior sample
    sEp_struct = spm_unvec(sEp(:,i),pE_struct);
    
    g = feval(PRF.M.IS, sEp_struct, PRF.M, PRF.U, 'get_response', xy);
    post_pd = post_pd + g;    
end

% Average
prior_pd = prior_pd .* (1/nsamp);
post_pd  = post_pd .* (1/nsamp);

spm_progress_bar('Clear');

% Close progress window
if ishandle(f_progress)
    close(f_progress);
end

% Return focus to main window
if ishandle(f)
    figure(f);
end