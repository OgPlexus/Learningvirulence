function period = preprocessCOVIDPeriod(Mu, Beta, startIdx, endIdx, periodId, periodLabel)
% preprocessCOVIDPeriod  Select one study period, drop non-finite pairs and min-max scale Beta within the period.
%
% CALLS     None
% CALLED BY analyzeCOVIDRegionLETR, buildSensitivitySystems
%
% INPUT
%   Mu, Beta       numeric column vectors of equal length (full record of one region)
%   startIdx, endIdx   first and last row of the period
%   periodId, periodLabel   char identifiers, for example 'preVacc' and 'Pre-mass-vaccination period'
% OUTPUT
%   period   struct with id, label, N, Mu, Beta (scaled to [0,1]), betaScaling (rawMin, rawMax),
%            ranges (mu and beta, each [min max])
MuCase = Mu(startIdx:endIdx);
BetaCase = Beta(startIdx:endIdx);
valid = isfinite(MuCase) & isfinite(BetaCase);
MuCase = MuCase(valid);
BetaCase = BetaCase(valid);
betaRawMin = min(BetaCase);
betaRawMax = max(BetaCase);
BetaScaled = (BetaCase - betaRawMin) ./ (betaRawMax - betaRawMin);
period.id = periodId;
period.label = periodLabel;
period.N = numel(MuCase);
period.Mu = MuCase;
period.Beta = BetaScaled;
period.betaScaling.rawMin = betaRawMin;
period.betaScaling.rawMax = betaRawMax;
period.ranges.mu = [min(MuCase), max(MuCase)];
period.ranges.beta = [min(BetaScaled), max(BetaScaled)];
end
