function cmap = batlow(N)
% BATLOW Generates the perceptually uniform, CVD-safe batlow colormap.
%   cmap = batlow(N) returns an N-by-3 matrix containing the RGB colormap.
%   Default N = 256.

if nargin < 1
    N = 256; 
end

% Anchor knot points for Fabio Crameri's batlow colormap (normalized RGB)
knots = [
    0.0118, 0.1176, 0.3804;  % Deep Navy
    0.0627, 0.2980, 0.5216;  % Blue
    0.0784, 0.4784, 0.5412;  % Teal
    0.2745, 0.6196, 0.4902;  % Muted Green
    0.6118, 0.7098, 0.3843;  % Light Olive Green
    0.9333, 0.7412, 0.3843;  % Warm Yellow
    0.9412, 0.5843, 0.5647;  % Salmon/Pink
    0.8667, 0.4157, 0.6588   % Pink-Magenta
];

x_knots = linspace(0, 1, size(knots, 1));
x_out   = linspace(0, 1, N)';

cmap = interp1(x_knots, knots, x_out, 'pchip');
end