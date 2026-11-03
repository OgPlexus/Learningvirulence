function saveCOVIDRegionResults(result, opts)
% saveCOVIDRegionResults  Save the result of one region as a MAT file and one GeoC table (CSV) per period.
%
% CALLS     None
% CALLED BY runCOVID
%
% INPUT
%   result   struct from analyzeCOVIDRegionLETR
%   opts     struct with OutputRoot (folder, created if missing for the region subfolder)
regionDir = fullfile(opts.OutputRoot, result.id);
if ~exist(regionDir,'dir'), mkdir(regionDir); end
save(fullfile(regionDir, sprintf('LETR_%s_Results.mat', result.id)), 'result', '-v7.3');
for p = 1:numel(result.periods)
    writetable(result.periods(p).GeoC.table, fullfile(regionDir, ...
        sprintf('LETR_%s_%s_GeoC.csv', result.id, result.periods(p).id)));
end
end
