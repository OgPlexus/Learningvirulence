function out = estimateInvariantDensity2D(mu, beta, varargin)
% estimateInvariantDensity2D  Trait transfer operator (boundary-reflected kernel density, sparse) and its invariant density.
%
% CALLS     None (uses sparse eigs from base MATLAB)
% CALLED BY runCase1a, runCase1b, runCase1c, analyzeCOVIDRegionLETR
%
% INPUT
%   mu, beta   numeric vectors of equal length, a time series of the two traits
%   Name-value options
%     NBins      grid cells per trait (default 100)
%     Bounds     [muMin betaMin; muMax betaMax], default is the data range with 2 percent padding
%     Bandwidth  'silverman' (default) or a two-element vector [h_mu h_beta]
%     Truncation kernel truncation radius in bandwidths (default 6)
%
% OUTPUT
%   out   struct with invariantDensity (NBins-by-NBins, rows index beta and columns index mu), MU and BETA
%         (meshgrid of cell centers), muVec, betaVec, f_mu_marginal, f_beta_marginal, bandwidth, bounds

p = inputParser;
addParameter(p, 'NBins', 100, @isnumeric);
addParameter(p, 'Bounds', [], @isnumeric);
addParameter(p, 'Bandwidth', 'silverman');
addParameter(p, 'Truncation', 6, @isnumeric);
parse(p, varargin{:});

nbins    = p.Results.NBins;
bounds   = p.Results.Bounds;
bwOpt    = p.Results.Bandwidth;
trunc    = p.Results.Truncation;

mu = mu(:); beta = beta(:);
valid = isfinite(mu) & isfinite(beta);
mu = mu(valid); beta = beta(valid);
mu_src = mu(1:end-1);   mu_tgt = mu(2:end);
beta_src = beta(1:end-1); beta_tgt = beta(2:end);

valid2 = isfinite(mu_src) & isfinite(mu_tgt) & isfinite(beta_src) & isfinite(beta_tgt);
mu_src = mu_src(valid2); mu_tgt = mu_tgt(valid2);
beta_src = beta_src(valid2); beta_tgt = beta_tgt(valid2);
N_T = numel(mu_src);

if isempty(bounds)
    padMu = 0.02 * (max(mu_src) - min(mu_src));
    padBeta = 0.02 * (max(beta_src) - min(beta_src));
    muMin = max(0, min(mu_src) - padMu);     muMax = max(mu_src) + padMu;
    betaMin = max(0, min(beta_src) - padBeta); betaMax = max(beta_src) + padBeta;
else
    muMin = bounds(1,1); betaMin = bounds(1,2);
    muMax = bounds(2,1); betaMax = bounds(2,2);
end

dMu   = (muMax - muMin) / nbins;
dBeta = (betaMax - betaMin) / nbins;
dA    = dMu * dBeta;

edges_mu   = linspace(muMin, muMax, nbins + 1);
edges_beta = linspace(betaMin, betaMax, nbins + 1);
muVec   = edges_mu(1:end-1) + dMu / 2;
betaVec = edges_beta(1:end-1) + dBeta / 2;
[MU, BETA] = meshgrid(muVec, betaVec);

if (ischar(bwOpt) || isstring(bwOpt)) && strcmpi(bwOpt, 'silverman')
    h_mu   = 1.06 * std(mu_tgt)   * N_T^(-1/5);
    h_beta = 1.06 * std(beta_tgt) * N_T^(-1/5);
else
    h_mu   = bwOpt(1);
    h_beta = bwOpt(2);
end

bin_mu   = min(max(discretizeLocal(mu_src, edges_mu), 1), nbins);
bin_beta = min(max(discretizeLocal(beta_src, edges_beta), 1), nbins);
cell_src = bin_beta + (bin_mu - 1) * nbins;

M = nbins * nbins;

rowIdxAll = cell(M, 1);
colIdxAll = cell(M, 1);
valAll    = cell(M, 1);
emptyRows = false(M, 1);

for i = 1:M
    S_i = find(cell_src == i);
    N_i = numel(S_i);
    if N_i == 0
        emptyRows(i) = true;
        continue;
    end

    mu_next   = mu_tgt(S_i);
    beta_next = beta_tgt(S_i);

    loMu = min([mu_next; 2*muMin - mu_next; 2*muMax - mu_next]) - trunc*h_mu;
    hiMu = max([mu_next; 2*muMin - mu_next; 2*muMax - mu_next]) + trunc*h_mu;
    keepMu = find(muVec >= loMu & muVec <= hiMu);
    if isempty(keepMu)
        [~, nearest] = min(abs(muVec - mean(mu_next)));
        keepMu = nearest;
    end

    loBeta = min([beta_next; 2*betaMin - beta_next; 2*betaMax - beta_next]) - trunc*h_beta;
    hiBeta = max([beta_next; 2*betaMin - beta_next; 2*betaMax - beta_next]) + trunc*h_beta;
    keepBeta = find(betaVec >= loBeta & betaVec <= hiBeta);
    if isempty(keepBeta)
        [~, nearest] = min(abs(betaVec - mean(beta_next)));
        keepBeta = nearest;
    end

    cMu = muVec(keepMu)';
    d_mu_1 = cMu - mu_next(:)';
    d_mu_2 = cMu - (2*muMin - mu_next(:)');
    d_mu_3 = cMu - (2*muMax - mu_next(:)');
    G_mu = (exp(-0.5*(d_mu_1/h_mu).^2) + exp(-0.5*(d_mu_2/h_mu).^2) + exp(-0.5*(d_mu_3/h_mu).^2)) ...
        / (sqrt(2*pi) * h_mu);
    G_mu = sum(G_mu, 2) / N_i;

    cBeta = betaVec(keepBeta)';
    d_beta_1 = cBeta - beta_next(:)';
    d_beta_2 = cBeta - (2*betaMin - beta_next(:)');
    d_beta_3 = cBeta - (2*betaMax - beta_next(:)');
    G_beta = (exp(-0.5*(d_beta_1/h_beta).^2) + exp(-0.5*(d_beta_2/h_beta).^2) + exp(-0.5*(d_beta_3/h_beta).^2)) ...
        / (sqrt(2*pi) * h_beta);
    G_beta = sum(G_beta, 2) / N_i;

    K_local = G_beta * G_mu';
    [KEEP_BETA, KEEP_MU] = ndgrid(keepBeta, keepMu);
    localCols = KEEP_BETA(:) + (KEEP_MU(:) - 1) * nbins;
    localVals = dA * K_local(:);
    s = sum(localVals);
    if s > 0
        localVals = localVals / s;
    else
        localVals(:) = 1 / numel(localVals);
    end

    rowIdxAll{i} = repmat(i, numel(localCols), 1);
    colIdxAll{i} = localCols;
    valAll{i}    = localVals;
end

rowIdx = cat(1, rowIdxAll{:});
colIdx = cat(1, colIdxAll{:});
vals   = cat(1, valAll{:});
if any(emptyRows)
    emptyList = find(emptyRows);
    nEmpty = numel(emptyList);
    rowIdx = [rowIdx; repelem(emptyList, M)];
    colIdx = [colIdx; repmat((1:M)', nEmpty, 1)];
    vals   = [vals;   repmat(1/M, nEmpty*M, 1)];
end
P = sparse(rowIdx, colIdx, vals, M, M);

pi_vec = solveStationarySparse(P);

pi_2D = reshape(pi_vec, [nbins, nbins]);
f2D   = pi_2D / dA;
f_mu_marginal   = sum(pi_2D, 1)' / dMu;
f_beta_marginal = sum(pi_2D, 2)   / dBeta;

out = struct();
out.invariantDensity = f2D;
out.f_mu_marginal    = f_mu_marginal;
out.f_beta_marginal  = f_beta_marginal;
out.MU               = MU;
out.BETA             = BETA;
out.muVec            = muVec;
out.betaVec          = betaVec;
out.bandwidth        = [h_mu, h_beta];
out.bounds           = [muMin, betaMin; muMax, betaMax];

end

function bin = discretizeLocal(x, edges)
n = numel(edges) - 1;
[~, ~, bin] = histcounts(x, edges);
bin(bin == 0) = NaN;
bin = min(bin, n);
end

function pi_vec = solveStationarySparse(P)
M = size(P, 1);
usedFallback = false;
try
    opts.tol = 1e-10;
    opts.maxit = 500;
    [V, D] = eigs(P', 1, 1, opts);
    pi_vec = real(V);
    if abs(D - 1) > 1e-3
        warning('estimateInvariantDensity2D:eigsEigenvalue', ...
            'eigs returned eigenvalue %.6f instead of 1; falling back to sparse power iteration.', real(D));
        usedFallback = true;
    end
catch ME
    warning('estimateInvariantDensity2D:eigsFailed', ...
        'eigs failed (%s); falling back to sparse power iteration.', ME.message);
    usedFallback = true;
end

if usedFallback
    pi_vec = ones(M, 1) / M;
    Pt = P';
    for iter = 1:5000
        pi_next = Pt * pi_vec;
        if norm(pi_next - pi_vec, 1) < 1e-12
            pi_vec = pi_next;
            break;
        end
        pi_vec = pi_next;
    end
end

pi_vec = max(0, pi_vec);
pi_vec = pi_vec / sum(pi_vec);
end
