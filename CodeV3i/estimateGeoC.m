function [T, details] = estimateGeoC(Mu, Beta, varargin)
% estimateGeoC  Geometric information flow (GeoC) between two series, in both directions.
%
% CALLS     conditionalD2
% CALLED BY runCase1a, runCase1b, runCase1c, analyzeCOVIDRegionLETR, runSampleSize
%
% INPUT
%   Mu, Beta   numeric vectors of equal length (virulence and transmission series)
%   Name-value options
%     Subsample  number of one-step transitions drawn at random (default 2500)
%     Seed       random seed for the subsample (default 1)
%     any other option is passed to estimateD2 (for example R2min, Cmax, SlopeCVmax, MaxSlopeDrift,
%     MinPlateauLength, MinPairs)
%
% OUTPUT
%   T          table with columns Term (names) and RawValue (D2 values and raw GeoC values),
%              including GeoC_beta_to_mu and GeoC_mu_to_beta
%   details    struct of estimateD2 outputs with fields D2_mun, D2_mun_munp1, D2_mun_betan,
%              D2_mun_betan_munp1, D2_betan, D2_betan_betanp1, D2_betan_mun, D2_betan_mun_betanp1
p = inputParser;
addParameter(p, 'Subsample', 2500);
addParameter(p, 'Seed', 1);
p.KeepUnmatched = true;
parse(p, varargin{:});
d2opts = namedargs2cell(p.Unmatched);

Mu = Mu(:); Beta = Beta(:);
xn = Mu(1:end-1); xnp1 = Mu(2:end);
yn = Beta(1:end-1); ynp1 = Beta(2:end);

rng(p.Results.Seed);
Nfull = length(xn);
Ns = min(p.Results.Subsample, Nfull);
idx = sort(randperm(Nfull, Ns));
xn = xn(idx); xnp1 = xnp1(idx); yn = yn(idx); ynp1 = ynp1(idx);

[~, Dcond_mu_raw, R_mu]            = conditionalD2(xn, xnp1, d2opts{:});
[~, Dcond_mu_beta_raw, R_mub] = conditionalD2([xn, yn], xnp1, d2opts{:});
GeoC_beta_to_mu = Dcond_mu_raw - Dcond_mu_beta_raw;

[~, Dcond_beta_raw, R_beta]        = conditionalD2(yn, ynp1, d2opts{:});
[~, Dcond_beta_mu_raw, R_betam] = conditionalD2([yn, xn], ynp1, d2opts{:});
GeoC_mu_to_beta = Dcond_beta_raw - Dcond_beta_mu_raw;

Term = {'D2(mu_n)'; 'D2([mu_n,mu_{n+1}])'; 'GeoC(mu_{n+1}|mu_n)'; ...
        'D2([mu_n,beta_n])'; 'D2([mu_n,beta_n,mu_{n+1}])'; 'GeoC(mu_{n+1}|mu_n,beta_n)'; ...
        'GeoC_beta_to_mu'; ...
        'D2(beta_n)'; 'D2([beta_n,beta_{n+1}])'; 'GeoC(beta_{n+1}|beta_n)'; ...
        'D2([beta_n,mu_n])'; 'D2([beta_n,mu_n,beta_{n+1}])'; 'GeoC(beta_{n+1}|beta_n,mu_n)'; ...
        'GeoC_mu_to_beta'};

RawValue = [R_mu.RY.D2; R_mu.RXY.D2; Dcond_mu_raw; ...
            R_mub.RY.D2; R_mub.RXY.D2; Dcond_mu_beta_raw; ...
            GeoC_beta_to_mu; ...
            R_beta.RY.D2; R_beta.RXY.D2; Dcond_beta_raw; ...
            R_betam.RY.D2; R_betam.RXY.D2; Dcond_beta_mu_raw; ...
            GeoC_mu_to_beta];

T = table(Term, RawValue);

details.D2_mun               = R_mu.RY;
details.D2_mun_munp1         = R_mu.RXY;
details.D2_mun_betan         = R_mub.RY;
details.D2_mun_betan_munp1   = R_mub.RXY;
details.D2_betan             = R_beta.RY;
details.D2_betan_betanp1     = R_beta.RXY;
details.D2_betan_mun         = R_betam.RY;
details.D2_betan_mun_betanp1 = R_betam.RXY;
end
