function R = transferEntropyKSG(src, tgt, varargin)
% transferEntropyKSG  One-step transfer entropy src to tgt (nats) with a k-nearest-neighbor estimator and a permutation test.
%
% CALLS     None (uses KDTreeSearcher, rangesearch from the Statistics and Machine Learning Toolbox)
% CALLED BY runCase1a
%
% INPUT
%   src, tgt   numeric vectors of equal length
%   Name-value options  Subsample (2000), K (4), NumSurrogates (20), Seed (1)
% OUTPUT
%   R   struct with TE, p (permutation p value), surrogateMean, K, N
p = inputParser;
addParameter(p,'Subsample',2000); addParameter(p,'K',4);
addParameter(p,'NumSurrogates',20); addParameter(p,'Seed',1);
parse(p,varargin{:});
K = p.Results.K; nS = p.Results.NumSurrogates;

src = src(:); tgt = tgt(:);
s = src(1:end-1); t0 = tgt(1:end-1); t1 = tgt(2:end);
rng(p.Results.Seed);
Ns  = min(p.Results.Subsample, numel(s));
idx = sort(randperm(numel(s), Ns));
s = s(idx); t0 = t0(idx); t1 = t1(idx);
zs = @(v) (v - mean(v)) ./ std(v);
s = zs(s); t0 = zs(t0); t1 = zs(t1);

R.TE = cmiKSG(t1, s, t0, K);
sur = zeros(nS,1);
for i = 1:nS
    sur(i) = cmiKSG(t1, s(randperm(Ns)), t0, K);
end
R.p = (1 + sum(sur >= R.TE)) / (nS + 1);
R.surrogateMean = mean(sur);
R.K = K; R.N = Ns;
end

function c = cmiKSG(x, y, z, k)
XYZ = [x y z];
mdl = KDTreeSearcher(XYZ, 'Distance', 'chebychev');
[~, d] = knnsearch(mdl, XYZ, 'K', k+1);
eps_i = d(:, end) * (1 - 1e-10);
nxz = countWithin([x z], eps_i);
nyz = countWithin([y z], eps_i);
nz  = countWithin(z,     eps_i);
c = mean(psi(k) - psi(nxz+1) - psi(nyz+1) + psi(nz+1));
end

function n = countWithin(A, eps_i)
mdl = KDTreeSearcher(A, 'Distance', 'chebychev');
n = zeros(size(A,1),1);
for i = 1:size(A,1)
    idx = rangesearch(mdl, A(i,:), eps_i(i));
    n(i) = numel(idx{1}) - 1;
end
end
