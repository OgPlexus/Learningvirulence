function [Dcond, Draw, R] = conditionalD2(Y, X, varargin)
% conditionalD2  Conditional dimension D2(X|Y) = D2([Y,X]) - D2(Y).
%
% CALLS     estimateD2
% CALLED BY estimateGeoC, estimateGeoCTI
%
% INPUT
%   Y, X       N-by-dy and N-by-dx numeric matrices with matching rows
%   varargin   name-value options passed to estimateD2
%
% OUTPUT
%   Dcond      max(0, Draw)
%   Draw       D2([Y,X]) - D2(Y), NaN when a plateau is missing
%   R          struct with RY, RXY (estimateD2 outputs) and success (both plateaus found)
RY  = estimateD2(Y, varargin{:});
RXY = estimateD2([Y, X], varargin{:});
Draw  = RXY.D2 - RY.D2;
Dcond = max(0, Draw);
R.RY = RY;
R.RXY = RXY;
R.success = RY.success && RXY.success;
if ~R.success
    warning('conditionalD2:noPlateau', 'No acceptable scaling region for D2(Y) or D2([Y,X]).');
end
end
