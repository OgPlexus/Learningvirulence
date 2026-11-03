function result = analyzeCOVIDRegionLETR(region, opts)
% analyzeCOVIDRegionLETR  GeoC and invariant density of one region for the full, pre- and post-mass-vaccination periods.
%
% CALLS     preprocessCOVIDPeriod, estimateGeoC, estimateInvariantDensity2D
% CALLED BY runCOVID
%
% INPUT
%   region   one element of defineCOVIDRegionConfigs (id, label, dataFile, periodBreakIndex)
%   opts     struct with
%              GeoC                struct whose fields are name-value options of estimateGeoC
%                                  (Subsample, Seed, MinPlateauLength, SlopeCVmax, MaxSlopeDrift, R2min, Cmax, MinPairs)
%              InvariantDensity    struct with NBins
% OUTPUT
%   result   struct with id, label, raw (Mu, Beta) and periods(1:3). Each period has id, label, N, Mu, Beta,
%            betaScaling, ranges, GeoC (table, details) and invariantDensity.
S = load(region.dataFile);
Mu = S.Mu(:);
Beta = S.Beta(:);
Ntotal = numel(Mu);
brk = region.periodBreakIndex;

periodDefs = { ...
    'full',     'Full period',                  1,   Ntotal; ...
    'preVacc',  'Pre-mass-vaccination period',  1,   brk-1; ...
    'postVacc', 'Post-mass-vaccination period', brk, Ntotal};

result.id = region.id;
result.label = region.label;
result.raw.Mu = Mu;
result.raw.Beta = Beta;

geoArgs = namedargs2cell(opts.GeoC);
for p = 1:size(periodDefs,1)
    period = preprocessCOVIDPeriod(Mu, Beta, periodDefs{p,3}, periodDefs{p,4}, periodDefs{p,1}, periodDefs{p,2});
    fprintf('%s, %s: N = %d\n', region.label, period.label, period.N);
    [period.GeoC.table, period.GeoC.details] = estimateGeoC(period.Mu, period.Beta, geoArgs{:});
    disp(period.GeoC.table);
    period.invariantDensity = estimateInvariantDensity2D(period.Mu, period.Beta, ...
        'NBins', opts.InvariantDensity.NBins, ...
        'Bounds', [min(period.Mu), min(period.Beta); max(period.Mu), max(period.Beta)]);
    result.periods(p) = period;
end
end
