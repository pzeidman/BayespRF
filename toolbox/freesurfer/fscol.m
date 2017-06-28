function cmap = fscol(res)
%
% cmap = fscol([res=360])
%
% Returns the FreeSurfer colour map (like MatLab colour maps) with the
% resolution res, which must be a multiple of 4. By default res is 360. 

if nargin == 0
    res = 360;   % steps from one colour to the next
end
steps = res / 4;

% Colour peaks
R = [1 0 0];
Y = [1 1 0];
G = [0 1 0];
B = [0 0 1];

% Create the colour map
cmap = [];
cmap = [cmap; linspace(R(1), Y(1), steps)', linspace(R(2), Y(2), steps)', linspace(R(3), Y(3), steps)'];
cmap = [cmap; linspace(Y(1), G(1), steps)', linspace(Y(2), G(2), steps)', linspace(Y(3), G(3), steps)'];
cmap = [cmap; linspace(G(1), B(1), steps)', linspace(G(2), B(2), steps)', linspace(G(3), B(3), steps)'];
cmap = [cmap; linspace(B(1), R(1), steps)', linspace(B(2), R(2), steps)', linspace(B(3), R(3), steps)'];
