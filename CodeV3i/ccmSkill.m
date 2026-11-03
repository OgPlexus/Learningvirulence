function R = ccmSkill(x, y, varargin)
% ccmSkill  Convergent cross mapping skill (simplex projection), x estimated from the shadow manifold of y.
%
% CALLS     None (uses KDTreeSearcher, knnsearch from the Statistics and Machine Learning Toolbox)
% CALLED BY runCase1a
%
% INPUT
%   x, y   numeric vectors of equal length
%   Name-value options  E (embedding dimension, 3), Tau (1), LibSizes ([100 300 1000 3000]),
%                       NumLibs (10), NumPred (300), Seed (1)
% OUTPUT
%   R   struct with L (library sizes used), rho, rhoSD, deltaRho (skill at largest minus smallest library), E, Tau
p = inputParser;
addParameter(p,'E',3); addParameter(p,'Tau',1);
addParameter(p,'LibSizes',[100 300 1000 3000]);
addParameter(p,'NumLibs',10); addParameter(p,'NumPred',300);
addParameter(p,'Seed',1);
parse(p,varargin{:});
E = p.Results.E; tau = p.Results.Tau; rng(p.Results.Seed);

x = x(:); y = y(:);
n  = numel(x) - (E-1)*tau;
My = zeros(n,E);
for i = 0:E-1
    My(:,i+1) = y((E-1-i)*tau + (1:n));
end
xs = x((E-1)*tau + (1:n));

Ls = p.Results.LibSizes(:)';
Ls = Ls(Ls < n - 1);
rho = nan(numel(Ls),1); rhoSD = nan(numel(Ls),1);
for a = 1:numel(Ls)
    L = Ls(a); r = nan(p.Results.NumLibs,1);
    for b = 1:p.Results.NumLibs
        i0   = randi(n - L);
        lib  = (i0:i0+L-1)';
        pred = randperm(n, min(p.Results.NumPred, n))';
        mdl  = KDTreeSearcher(My(lib,:));
        [ii, dd] = knnsearch(mdl, My(pred,:), 'K', E+2);
        libIdx = lib(ii);
        xhat = zeros(numel(pred),1);
        for q = 1:numel(pred)
            keep = find(libIdx(q,:) ~= pred(q), E+1);
            d = dd(q,keep);
            w = exp(-d ./ max(d(1), 1e-12));
            xhat(q) = sum(w .* xs(libIdx(q,keep))') / sum(w);
        end
        r(b) = corr(xhat, xs(pred));
    end
    rho(a) = mean(r, 'omitnan'); rhoSD(a) = std(r, 'omitnan');
end
R.L = Ls; R.rho = rho; R.rhoSD = rhoSD;
R.deltaRho = rho(end) - rho(1);
R.E = E; R.Tau = tau;
end
