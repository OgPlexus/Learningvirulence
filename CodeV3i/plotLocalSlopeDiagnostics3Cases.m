function fig = plotLocalSlopeDiagnostics3Cases(detailsCases, caseNames)
% plotLocalSlopeDiagnostics3Cases  3-by-4 local-slope figure, one row per case and one column per GeoC term.
%
% CALLS     drawLocalSlopePanel
% CALLED BY makeCOVIDRegionFigures
%
% INPUT
%   detailsCases   1x3 cell array of estimateGeoC details structs (rows of the figure)
%   caseNames      1x3 cell array of row labels
% OUTPUT
%   fig            figure handle
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 24, 17]);
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
columnTitles = { ...
    '$\mathrm{GeoC}(\mu_{n+1}\mid\mu_n)$', ...
    '$\mathrm{GeoC}(\mu_{n+1}\mid\mu_n,\beta_n)$', ...
    '$\mathrm{GeoC}(\beta_{n+1}\mid\beta_n)$', ...
    '$\mathrm{GeoC}(\beta_{n+1}\mid\beta_n,\mu_n)$'};
letters = 'abcdefghijkl';

for c = 1:3
    d = detailsCases{c};
    panels = {
        d.D2_mun,       d.D2_mun_munp1,        'D_2(\mu_n)',           'D_2([\mu_n,\mu_{n+1}])';
        d.D2_mun_betan, d.D2_mun_betan_munp1,  'D_2([\mu_n,\beta_n])', 'D_2([\mu_n,\beta_n,\mu_{n+1}])';
        d.D2_betan,     d.D2_betan_betanp1,    'D_2(\beta_n)',         'D_2([\beta_n,\beta_{n+1}])';
        d.D2_betan_mun, d.D2_betan_mun_betanp1,'D_2([\beta_n,\mu_n])', 'D_2([\beta_n,\mu_n,\beta_{n+1}])'};
    for k = 1:4
        nexttile;
        drawLocalSlopePanel(panels{k,1}, panels{k,2}, panels{k,3}, panels{k,4}, []);
        if k == 1, ylabel('$d\log C / d\log r$', 'Interpreter', 'latex', 'FontSize', 9); end
        if c == 3, xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 9); end
        if c == 1, title(columnTitles{k}, 'Interpreter', 'latex', 'FontSize', 9.5, 'FontWeight', 'normal'); end
        text(-0.12, 1.04, sprintf('\\textbf{(%s)}', letters((c - 1) * 4 + k)), ...
            'Units', 'normalized', 'FontSize', 9, 'Interpreter', 'latex');
        hold off;
    end
end

for c = 1:3
    nexttile((c - 1) * 4 + 4);
    text(1.08, 0.5, caseNames{c}, 'Units', 'normalized', 'Rotation', 270, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'Helvetica', 'FontSize', 9, 'Clipping', 'off');
end
set(fig, 'InvertHardcopy', 'off');
end
