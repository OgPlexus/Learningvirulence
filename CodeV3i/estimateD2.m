function R = estimateD2(X, varargin)
% estimateD2  Correlation dimension from a plateau of the local-slope profile of log C(r) against log r.
%
% CALLS     None (uses KDTreeSearcher, knnsearch, rangesearch from the Statistics and Machine Learning Toolbox)
% CALLED BY conditionalD2
%
% INPUT
%   X        N-by-d numeric matrix, one observation per row. Rows with non-finite values are dropped.
%   Name-value options (all optional)
%     NumRadii 80, MinPairs 1000, Cmax 0.20, MinPlateauLength 10, SlopeCVmax 0.15,
%     MaxSlopeDrift 0.10, R2min 0.98, LocalSlopeHalfWindow 2, TheilerWindow 0,
%     TimeIndex [] (time index of each row of X, used only when TheilerWindow > 0), RadiusK 10
%
% OUTPUT
%   R  struct with D2, R2, i, j (first and last selected radius index), r, C, logR, logC, localSlope,
%      localR2, slopeCV, slopeDrift, Npairs, success (true when a plateau was found).
%      D2 is NaN and success is false when no plateau is found.
opt = parseD2Options(varargin{:});
Noriginal = size(X,1);
finiteRows = all(isfinite(X), 2);
X = X(finiteRows, :);
N = size(X,1);
if N < 30
    R = emptyD2Result(); return;
end

if isempty(opt.TimeIndex)
    TimeIndex = (1:N)';
else
    TimeIndex = opt.TimeIndex(:);
    if numel(TimeIndex) == Noriginal
        TimeIndex = TimeIndex(finiteRows);
    end
end

for k = 1:size(X,2)
    s = std(X(:,k));
    if ~isfinite(s) || s <= eps
        R = emptyD2Result(); return;
    end
    X(:,k) = (X(:,k) - mean(X(:,k))) / s;
end

tree = KDTreeSearcher(X);
if opt.TheilerWindow > 0
    dt = abs(TimeIndex - TimeIndex');
    Npairs = sum(dt(triu(true(N), 1)) > opt.TheilerWindow);
else
    Npairs = N*(N-1)/2;
end

if Npairs < opt.MinPairs
    R = emptyD2Result(); R.Npairs = Npairs; R.TheilerWindow = opt.TheilerWindow; return;
end

K = min(opt.RadiusK+1, N);
[~, knnDist] = knnsearch(tree, X, 'K', K);
kDist = knnDist(:, K);
kDist = kDist(isfinite(kDist) & kDist > 0);
if isempty(kDist)
    R = emptyD2Result(); R.Npairs = Npairs; R.TheilerWindow = opt.TheilerWindow; return;
end

rmin = prctile(kDist, 5);
rmax = prctile(kDist, 95);
if ~isfinite(rmin) || ~isfinite(rmax) || rmax <= rmin
    R = emptyD2Result(); R.Npairs = Npairs; R.TheilerWindow = opt.TheilerWindow; return;
end

radii = logspace(log10(rmin), log10(rmax), opt.NumRadii);

pairCounts = zeros(size(radii));
for i = 1:N-1
    [idxCell, distCell] = rangesearch(tree, X(i,:), rmax);
    idx = idxCell{1}; dist = distCell{1};
    keep = idx > i;
    if opt.TheilerWindow > 0
        keep = keep & abs(TimeIndex(idx) - TimeIndex(i)) > opt.TheilerWindow;
    end
    dist = dist(keep);
    if isempty(dist), continue; end
    counts = histcounts(dist, [-Inf, radii]);
    pairCounts = pairCounts + cumsum(counts);
end
C = pairCounts / Npairs;

Cmin = max(opt.MinPairs/Npairs, 1/Npairs);
valid = C >= Cmin & C <= opt.Cmax & C > 0 & C < 1 & isfinite(C);
r = radii(valid); c = C(valid);
if numel(r) < opt.MinPlateauLength
    R = emptyD2Result(); R.radii = radii; R.C = C; R.r = r; R.c = c;
    R.Npairs = Npairs; R.Cmin = Cmin; R.Cmax = opt.Cmax; R.TheilerWindow = opt.TheilerWindow; return;
end

logR = log(r); logC = log(c); n = numel(logR);
h = opt.LocalSlopeHalfWindow;
localSlope = nan(n,1); localR2 = nan(n,1);

for k = 1:n
    i1 = max(1, k-h); i2 = min(n, k+h);
    if i2-i1+1 < 3, continue; end
    xx = logR(i1:i2); yy = logC(i1:i2);
    pp = polyfit(xx, yy, 1);
    yhat = polyval(pp, xx);
    ssres = sum((yy - yhat).^2);
    sstot = sum((yy - mean(yy)).^2);
    localSlope(k) = pp(1);
    if sstot > 0, localR2(k) = 1 - ssres/sstot; end
end

[candidateI, candidateJ] = findCandidatePlateaus(localSlope, logR, logC, opt);
if isempty(candidateI)
    R = emptyD2Result(); R.radii = radii; R.C = C; R.r = r; R.c = c;
    R.logR = logR; R.logC = logC; R.localSlope = localSlope; R.localR2 = localR2;
    R.Npairs = Npairs; R.Cmin = Cmin; R.Cmax = opt.Cmax; R.TheilerWindow = opt.TheilerWindow; return;
end

nCandidates = numel(candidateI);
expandedI = nan(nCandidates,1); expandedJ = nan(nCandidates,1);
candidateD2 = nan(nCandidates,1); candidateR2 = nan(nCandidates,1);
candidateCV = nan(nCandidates,1); candidateDrift = nan(nCandidates,1);
candidateLength = nan(nCandidates,1);

for q = 1:nCandidates
    [ii, jj] = expandPlateau(candidateI(q), candidateJ(q), localSlope, opt);
    [d2, r2, cv, drift] = fitWindow(logR, logC, localSlope, ii, jj);
    if r2 >= opt.R2min && cv <= opt.SlopeCVmax && drift <= opt.MaxSlopeDrift
        finalI = ii; finalJ = jj;
    else
        finalI = candidateI(q); finalJ = candidateJ(q);
        [d2, r2, cv, drift] = fitWindow(logR, logC, localSlope, finalI, finalJ);
    end
    expandedI(q) = finalI; expandedJ(q) = finalJ;
    candidateD2(q) = d2; candidateR2(q) = r2; candidateCV(q) = cv;
    candidateDrift(q) = drift; candidateLength(q) = finalJ - finalI + 1;
end

good = isfinite(candidateD2) & isfinite(candidateR2) & isfinite(candidateCV) & isfinite(candidateDrift);
if ~any(good)
    R = emptyD2Result(); R.radii = radii; R.C = C; R.r = r; R.c = c;
    R.logR = logR; R.logC = logC; R.localSlope = localSlope; R.localR2 = localR2;
    R.Npairs = Npairs; R.Cmin = Cmin; R.Cmax = opt.Cmax; R.TheilerWindow = opt.TheilerWindow; return;
end

goodIndices = find(good);
score = nan(numel(goodIndices), 1);
for q = 1:numel(goodIndices)
    k = goodIndices(q);
    score(q) = candidateLength(k) * candidateR2(k)^2 / (1 + candidateCV(k));
end
[~, bestLocal] = max(score);
best = goodIndices(bestLocal);

i = expandedI(best); j = expandedJ(best);
R.D2 = candidateD2(best); R.R2 = candidateR2(best);
R.i = i; R.j = j; R.windowLength = j - i + 1;
R.rMinSelected = r(i); R.rMaxSelected = r(j);
R.slopeCV = candidateCV(best); R.slopeDrift = candidateDrift(best);
R.Npairs = Npairs; R.Cmin = Cmin; R.Cmax = opt.Cmax; R.TheilerWindow = opt.TheilerWindow;
R.radii = radii; R.C = C; R.r = r; R.c = c; R.logR = logR; R.logC = logC;
R.localSlope = localSlope; R.localR2 = localR2; R.seedI = candidateI(best); R.seedJ = candidateJ(best);
R.success = true;
end

function opt = parseD2Options(varargin)
p = inputParser;
addParameter(p,'NumRadii',80); addParameter(p,'MinPairs',1000);
addParameter(p,'Cmax',0.20); addParameter(p,'MinPlateauLength',10);
addParameter(p,'SlopeCVmax',0.15); addParameter(p,'MaxSlopeDrift',0.10);
addParameter(p,'R2min',0.98); addParameter(p,'LocalSlopeHalfWindow',2);
addParameter(p,'TheilerWindow',0); addParameter(p,'TimeIndex',[]);
addParameter(p,'RadiusK',10);
parse(p, varargin{:}); opt = p.Results;
end

function [I, J] = findCandidatePlateaus(localSlope, logR, logC, opt)
n = numel(localSlope); I = []; J = [];
for a = 1:(n - opt.MinPlateauLength + 1)
    b = a + opt.MinPlateauLength - 1;
    ls = localSlope(a:b);
    if any(~isfinite(ls)) || any(ls <= 0), continue; end
    m = median(ls);
    if ~isfinite(m) || m <= 0, continue; end
    cv = std(ls)/m; drift = max(abs(ls-m))/m;
    if cv > opt.SlopeCVmax || drift > opt.MaxSlopeDrift, continue; end
    [~, R2, ~, ~] = fitWindow(logR, logC, localSlope, a, b);
    if ~isfinite(R2) || R2 < opt.R2min, continue; end
    I(end+1,1) = a; J(end+1,1) = b;
end
end

function [i, j] = expandPlateau(seedI, seedJ, localSlope, opt)
i = seedI; j = seedJ;
refSlope = median(localSlope(seedI:seedJ)); n = numel(localSlope);
while i > 1
    testSlope = localSlope(i-1:j);
    if any(~isfinite(testSlope)) || any(testSlope <= 0), break; end
    m = median(testSlope); if m <= 0, break; end
    if std(testSlope)/m <= opt.SlopeCVmax && max(abs(testSlope-refSlope))/max(abs(refSlope),eps) <= opt.MaxSlopeDrift
        i = i - 1;
    else, break; end
end
while j < n
    testSlope = localSlope(i:j+1);
    if any(~isfinite(testSlope)) || any(testSlope <= 0), break; end
    m = median(testSlope); if m <= 0, break; end
    if std(testSlope)/m <= opt.SlopeCVmax && max(abs(testSlope-refSlope))/max(abs(refSlope),eps) <= opt.MaxSlopeDrift
        j = j + 1;
    else, break; end
end
end

function [D2, R2, slopeCV, slopeDrift] = fitWindow(logR, logC, localSlope, i, j)
xx = logR(i:j); yy = logC(i:j);
if numel(xx) < 2 || any(~isfinite(xx)) || any(~isfinite(yy))
    D2 = NaN; R2 = NaN; slopeCV = Inf; slopeDrift = Inf; return;
end
pp = polyfit(xx, yy, 1);
yyhat = polyval(pp, xx);
SSres = sum((yy - yyhat).^2); SStot = sum((yy - mean(yy)).^2);
R2 = 1 - SSres / max(SStot, eps); D2 = pp(1);
ls = localSlope(i:j); m = median(ls);
if ~isfinite(m) || abs(m) <= eps
    slopeCV = Inf; slopeDrift = Inf; return;
end
slopeCV = std(ls) / abs(m); slopeDrift = max(abs(ls - m)) / abs(m);
end

function R = emptyD2Result()
R.D2 = NaN; R.R2 = NaN; R.i = NaN; R.j = NaN; R.windowLength = NaN;
R.rMinSelected = NaN; R.rMaxSelected = NaN; R.slopeCV = NaN; R.slopeDrift = NaN;
R.Npairs = NaN; R.Cmin = NaN; R.Cmax = NaN; R.TheilerWindow = NaN;
R.radii = []; R.C = []; R.r = []; R.c = []; R.logR = []; R.logC = [];
R.localSlope = []; R.localR2 = []; R.seedI = NaN; R.seedJ = NaN; R.success = false;
end
