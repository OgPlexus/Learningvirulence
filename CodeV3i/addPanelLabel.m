function addPanelLabel(ax, labelText, dx, dy, fs)
% addPanelLabel  Bold panel label such as (a) in normalized axes coordinates.
%
% CALLS     None
% CALLED BY runCase1a, runCase1b, runCase1c, plotInvariantDensity3Periods
%
% INPUT
%   ax         axes handle
%   labelText  char, for example '(a)'
%   dx, dy     (optional) normalized position, defaults -0.23 and 0.95
%   fs         (optional) font size, default 10.5
if nargin < 3, dx = -0.23; end
if nargin < 4, dy = 0.95; end
if nargin < 5, fs = 10.5; end
text(ax, dx, dy, ['\textbf{', labelText, '}'], 'Units', 'normalized', ...
    'Interpreter', 'latex', 'FontSize', fs, 'VerticalAlignment', 'bottom');
end
