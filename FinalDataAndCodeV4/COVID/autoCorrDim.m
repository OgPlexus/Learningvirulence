function [dim, params] = autoCorrDim(X, varargin)
% autoCorrDim  Simple automatic wrapper for MATLAB correlationDimension.
%   [dim, params] = autoCorrDim(X) computes an automatic NumPoints/MinRadius/
%   MaxRadius and returns the estimated correlation dimension 'dim' and a
%   params struct with fields .NumPoints, .MinRadius, .MaxRadius.
%
%   Optional name-value arguments:
%     'NumPoints'  - integer (default 100)
%     'MinRadius'  - numeric (if given, override auto-selection)
%     'MaxRadius'  - numeric (if given, override auto-selection)
%
% Example:
%   [d,p] = autoCorrDim(X1);
%   [d2,p2] = autoCorrDim([X1 X2], 'NumPoints',200);

% ---- parse inputs ----
p = inputParser;
addParameter(p, 'NumPoints', 100, @(x)isnumeric(x) && isscalar(x));
addParameter(p, 'MinRadius', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'MaxRadius', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});
Np = p.Results.NumPoints;
minR_user = p.Results.MinRadius;
maxR_user = p.Results.MaxRadius;

% ---- prepare data ----
if isvector(X)
    X = X(:); % column vector
end

% ---- get coarse correlation integral ----
% we request rRange and corInt so we can auto-select radii
try
    [~, rRange, corInt] = correlationDimension(X, [], 1, 'NumPoints', Np);
catch ME
    error('Error calling correlationDimension: %s', ME.message)
end

% clean nonpositive values (log needs positive)
valid = (rRange > 0) & (corInt > 0);
rRange = rRange(valid);
corInt = corInt(valid);

% ---- automatic min/max radius selection (unless user gave them) ----
if ~isempty(minR_user) && ~isempty(maxR_user)
    minR = minR_user;
    maxR = maxR_user;
else
    % fallback if too few points
    if numel(rRange) < 8
        % use empirical quantiles
        i1 = max(1, round(0.1*numel(rRange)));
        i2 = min(numel(rRange), round(0.9*numel(rRange)));
        minR = rRange(i1);
        maxR = rRange(i2);
    else
        lr = log(rRange(:));
        lc = log(corInt(:));
        win = min(6, floor(numel(lr)/2)); % sliding window size
        if win < 3
            i1 = max(1, round(0.1*numel(rRange)));
            i2 = min(numel(rRange), round(0.9*numel(rRange)));
            minR = rRange(i1);
            maxR = rRange(i2);
        else
            slopes = nan(numel(lr)-win+1,1);
            for i = 1:numel(slopes)
                pf = polyfit(lr(i:i+win-1), lc(i:i+win-1), 1);
                slopes(i) = pf(1);
            end
            % look for contiguous region where slope changes little
            ds = abs(diff(slopes));
            % relative threshold (tunable)
            thresh = 0.03 * max(1, median(abs(slopes)));
            goodIdx = find(ds < thresh);
            if isempty(goodIdx)
                % no plateau found: pick the window with most stable smoothed slope
                ssm = smooth(slopes,3);
                [~, mid] = max(ssm); % choose strongest slope region (heuristic)
                iStart = max(1, mid);
                iEnd = min(numel(lr), mid + win - 1);
            else
                iStart = goodIdx(1);
                iEnd = min(numel(lr), goodIdx(end) + win - 1);
            end
            minR = rRange(iStart);
            maxR = rRange(iEnd);
        end
    end
    % if user gave only one of min/max, allow override
    if ~isempty(minR_user), minR = minR_user; end
    if ~isempty(maxR_user), maxR = maxR_user; end
end

% safety checks: ensure minR < maxR and values inside observed range
if minR <= 0 || maxR <= 0 || minR >= maxR
    % choose mid-range fallback
    minR = rRange(max(1, floor(0.1*numel(rRange))));
    maxR = rRange(min(numel(rRange), floor(0.9*numel(rRange))));
end

% ---- recompute correlation dimension using chosen radii ----
dim = correlationDimension(X, [], 1, 'NumPoints', Np, 'MinRadius', minR, 'MaxRadius', maxR);

% ---- output params ----
params = struct('NumPoints', Np, 'MinRadius', minR, 'MaxRadius', maxR);

end
