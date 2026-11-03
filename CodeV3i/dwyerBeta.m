function beta = dwyerBeta(mu,delta,q0,q1)
% dwyerBeta  Transmission of a strain with virulence mu from the within-host titer model of Dwyer et al.
%
% CALLS     None
% CALLED BY generateCase1Series, runCase1c
%
% INPUT
%   mu              scalar virulence in (0,1)
%   delta, q0, q1   scalars, transmission scale and the two coefficients of the infectivity map
% OUTPUT
%   beta            scalar transmission
mu = min(max(mu,1e-12),1-1e-12);
k1 = 3.68 - 7.8.*mu + 7.5.*mu.^2;
k2 = max(0.265 - 0.66.*mu + 0.56.*mu.^2, 1e-12);
T  = max(334.7 - 70.3.*log(100.*mu), 0);
vfun = @(t) k1.*t.*exp(-k2.*t);
qfun = @(t) min(1,max(0,q0 + q1.*vfun(t)));
beta = max(delta.*integral(qfun,0,T,'ArrayValued',true), 0);
end
