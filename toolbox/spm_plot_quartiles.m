function [h,h2] = spm_plot_quartiles(E,x,varargin)
% Plot median (as a circle) and inter-quartile range (as a bar)
%
% E       - data matrix [observations x data points]
% x       - domain (x-axis)
% options - Plot options
%
% Returns:
%
% h       - handles of the data points
% h2      - handles of the error bars
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

if nargin < 2 || isempty(x)
    x = 1:size(E,2);
end

h  = [];
h2 = [];

for k = 1:size(E,2)
    
    data = E(:,k);
    
    % Compute statistics 
    e_mid = median(data);
    n     = length(data);

    % Calculate Tukey's hinges
    if mod(n,2) == 0
        e1 = data(data < e_mid);
        e2 = data(data > e_mid);
    else
        e1 = data(data <= e_mid);
        e2 = data(data >= e_mid);    
    end
    q1 = median(e1);
    q3 = median(e2);
    
    % Plot lower and upper whiskers
    h2(end+1) = line([x(k) x(k)],[q1 q3],...
        'LineWidth',1,'Color','k','LineWidth',2,varargin{:});
    hold on;
    
    % Plot median
    h(end+1) = plot(x(k),e_mid,'.','MarkerSize',20,'Color','k',varargin{:});        
end
hold off;