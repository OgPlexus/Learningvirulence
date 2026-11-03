function out = estimateGeoCTI(Mu, Beta, varargin)
% estimateGeoCTI  Compact GeoC estimate that returns numbers only and passes time indices to estimateD2.
%
% CALLS     conditionalD2
% CALLED BY runSensitivity, runSampleSize
%
% INPUT
%   Mu, Beta   numeric vectors of equal length
%   Name-value options
%     Subsample (default 2500), Seed (default 1), and any estimateD2 option such as TheilerWindow,
%     R2min, SlopeCVmax, MaxSlopeDrift, Cmax, NumRadii, LocalSlopeHalfWindow, RadiusK
%
% OUTPUT
%   out   struct with G_mu, G_mu_beta, G_beta, G_beta_mu (raw conditional GeoC values), B2M (beta to mu),
%         M2B (mu to beta), nFail (number of D2 pairs without a plateau), N (transitions used)
p = inputParser; p.KeepUnmatched = true;
addParameter(p,'Subsample',2500); addParameter(p,'Seed',1);
parse(p,varargin{:});
d2opts = namedargs2cell(p.Unmatched);

Mu = Mu(:); Beta = Beta(:);
xn = Mu(1:end-1);   xp = Mu(2:end);
yn = Beta(1:end-1); yp = Beta(2:end);
rng(p.Results.Seed);
Nfull = numel(xn);
Ns  = min(p.Results.Subsample, Nfull);
idx = sort(randperm(Nfull, Ns));
xn = xn(idx); xp = xp(idx); yn = yn(idx); yp = yp(idx);
d2opts = [d2opts, {'TimeIndex', idx(:)}];

ws = warning('off','conditionalD2:noPlateau'); cleanupObj = onCleanup(@() warning(ws));
[~, r1, R1] = conditionalD2(xn,       xp, d2opts{:});
[~, r2, R2] = conditionalD2([xn yn],  xp, d2opts{:});
[~, r3, R3] = conditionalD2(yn,       yp, d2opts{:});
[~, r4, R4] = conditionalD2([yn xn],  yp, d2opts{:});

out.G_mu = r1; out.G_mu_beta = r2; out.G_beta = r3; out.G_beta_mu = r4;
out.B2M = r1 - r2;
out.M2B = r3 - r4;
out.nFail = 4 - (R1.success + R2.success + R3.success + R4.success);
out.N = Ns;
end
