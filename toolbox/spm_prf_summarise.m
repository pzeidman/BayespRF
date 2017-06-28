function y = spm_prf_summarise(PRF,xY,method,noplot)
% Summarises and optionally plots a pRF receptive field within an ROI
%
% PRF    - PRF file or mat
% xY     - mask image
% method - 'sum' or 'BPA' (default: sum)
% noplot - if true, disables plotting (default: false)
% 
% y      - matrix representation of the receptive field
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

if ischar(PRF)
    PRF = load(PRF);
    PRF = PRF.PRF;
end

if nargin < 3 || isempty(method)
    method = 'sum';
end

if nargin < 4
    noplot = false;
end

% List X,Y stimulus coordinates
pmax = PRF.M.pmax/2;
x_bins = -pmax:0.5:pmax;
y_bins = x_bins;
b      = length(x_bins);
[x2,y2] = meshgrid(x_bins,y_bins);
xy = [x2(:) y2(:)];

 % ROI definition
if ~isstruct(xY)
    xY = struct('def','mask','spec',xY);
end

if strcmpi(method,'sum')
    % Summed prF responses
    if size(PRF.Y.y,2) == 1
        included_voxels = 1;
    else
        [~, ~, included_voxels] = spm_ROI(xY, PRF.xY.XYZmm);
    end
    roi_sum = [];
    for i = 1:length(included_voxels)        
        g = feval(PRF.M.IS, PRF.Ep{included_voxels(i)}, PRF.M, PRF.U, 'get_response', xy);

        if i == 1
            roi_sum = g;
        else
            roi_sum = roi_sum + g;
        end
    end   
    y = roi_sum;
elseif strcmpi(method,'BPA')
    % BPA over voxels
    BPA = spm_prf_bpa_within_subject(PRF,xY,true);
    y = feval(PRF.M.IS, BPA.Ep{1}, PRF.M, PRF.U, 'get_response', xy);
else
    error('Unknown summarise method');
end

% Reshape to 2x2 stim space
y = reshape(y,b,b);  

if noplot, return; end

% Plot
mid = ceil(b/2);
ticks = [1 mid b];   
imagesc(y); axis square;
set(gca,'YDir','normal');
set(gca,'XTickLabel',x_bins(ticks),'YTickLabel',x_bins(ticks),'XTick',ticks,'YTick',ticks); 
colormap(jet()); hold on;    
line([1 b],[mid mid],'Color','w');
line([mid mid],[1 b],'Color','w');
drawnow(); 