% runCOVID  LETR analysis of SARS-CoV-2 data for all regions (full, pre- and post-mass-vaccination periods).
%
% CALLS     defineCOVIDRegionConfigs, analyzeCOVIDRegionLETR, saveCOVIDRegionResults, makeCOVIDRegionFigures,
%           makeCOVIDRegionalMainFigure
% INPUT     MAT files listed by defineCOVIDRegionConfigs in the COVID folder, each with column vectors Mu and Beta
% OUTPUT    ResultsCOVIDMultiRegion/<region>/ (MAT file, GeoC tables, two figures per region),
%           ResultsCOVIDMultiRegion/MainFigures/FigureRegionalCOVIDMain.png and LETRCOVIDAllRegionsResults.mat

clear; close all; clc;
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

opts.OutputRoot = 'ResultsCOVIDMultiRegion';
opts.GeoC = struct('Subsample', 2500, 'Seed', 2, 'MinPlateauLength', 10, 'SlopeCVmax', 0.20, ...
    'MaxSlopeDrift', 0.15, 'R2min', 0.95, 'Cmax', 0.20, 'MinPairs', 1000);
opts.InvariantDensity.NBins = 100;
opts.Figures.Resolution = 600;
opts.Figures.ExportPDF = true;

%%
regions = defineCOVIDRegionConfigs();
if ~exist(opts.OutputRoot, 'dir'), mkdir(opts.OutputRoot); end

results = cell(numel(regions), 1);
for r = 1:numel(regions)
    results{r} = analyzeCOVIDRegionLETR(regions(r), opts);
    saveCOVIDRegionResults(results{r}, opts);
    makeCOVIDRegionFigures(results{r}, opts);
end
results = vertcat(results{:});

%%
makeCOVIDRegionalMainFigure(results, opts);
save(fullfile(opts.OutputRoot, 'LETRCOVIDAllRegionsResults.mat'), 'results', 'regions', 'opts', '-v7.3');
