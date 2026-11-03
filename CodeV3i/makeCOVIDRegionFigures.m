function makeCOVIDRegionFigures(result, opts)
% makeCOVIDRegionFigures  Save the local-slope figure and the invariant-density figure of one region.
%
% CALLS     plotLocalSlopeDiagnostics3Cases, plotInvariantDensity3Periods
% CALLED BY runCOVID
%
% INPUT
%   result   struct from analyzeCOVIDRegionLETR
%   opts     struct with OutputRoot, Figures.Resolution (dpi) and Figures.ExportPDF (logical).
%            The folder OutputRoot/result.id must already exist (saveCOVIDRegionResults creates it).
regionDir = fullfile(opts.OutputRoot, result.id);
detailsCases = arrayfun(@(p) p.GeoC.details, result.periods, 'UniformOutput', false);
caseNames = {result.periods.label};

fig = plotLocalSlopeDiagnostics3Cases(detailsCases, caseNames);
exportgraphics(fig, fullfile(regionDir, ['Figure_' result.id '_GeoC_LocalSlopeDiagnostics.png']), ...
    'Resolution', opts.Figures.Resolution, 'BackgroundColor', 'white');
if opts.Figures.ExportPDF
    exportgraphics(fig, fullfile(regionDir, ['Figure_' result.id '_GeoC_LocalSlopeDiagnostics.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
end
close(fig);

fig = plotInvariantDensity3Periods(result);
exportgraphics(fig, fullfile(regionDir, ['Figure_' result.id '_InvariantDensity_ThreePeriods.png']), ...
    'Resolution', opts.Figures.Resolution, 'BackgroundColor', 'white');
if opts.Figures.ExportPDF
    exportgraphics(fig, fullfile(regionDir, ['Figure_' result.id '_InvariantDensity_ThreePeriods.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
end
close(fig);
end
