function fig = plotLocalSlopeDiagnostics2x2(details)
% plotLocalSlopeDiagnostics2x2  2-by-2 local-slope figure for the four GeoC terms of one record.
%
% CALLS     None
% CALLED BY runCase1a, runCase1b
%
% INPUT
%   details   estimateGeoC details struct
% OUTPUT
%   fig       figure handle

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 18, 13]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

c1 = [0.12, 0.31, 0.89];
c2 = [0.85, 0.15, 0.20];

panels = {
    details.D2_mun, details.D2_mun_munp1, ...
    '$D_2(\mu_n)$', '$D_2([\mu_n,\mu_{n+1}])$', ...
    '$\mathrm{GeoC}(\mu_{n+1}\mid\mu_n)$', 'a'; ...

    details.D2_betan, details.D2_betan_betanp1, ...
    '$D_2(\beta_n)$', '$D_2([\beta_n,\beta_{n+1}])$', ...
    '$\mathrm{GeoC}(\beta_{n+1}\mid\beta_n)$', 'b'; ...

    details.D2_mun_betan, details.D2_mun_betan_munp1, ...
    '$D_2([\mu_n,\beta_n])$', '$D_2([\mu_n,\beta_n,\mu_{n+1}])$', ...
    '$\mathrm{GeoC}(\mu_{n+1}\mid\mu_n,\beta_n)$', 'c'; ...

    details.D2_betan_mun, details.D2_betan_mun_betanp1, ...
    '$D_2([\beta_n,\mu_n])$', '$D_2([\beta_n,\mu_n,\beta_{n+1}])$', ...
    '$\mathrm{GeoC}(\beta_{n+1}\mid\beta_n,\mu_n)$', 'd'
};

for k = 1:4
    nexttile;
    hold on;

    R1   = panels{k, 1};
    R2   = panels{k, 2};
    lbl1 = panels{k, 3};
    lbl2 = panels{k, 4};
    pTag = panels{k, 6};

    hLeg = [];
    tLeg = {};

    plateauSlopes = [];
    if isstruct(R1) && R1.success
        idx1 = max(1, R1.i-3) : min(length(R1.localSlope), R1.j+3);
        plateauSlopes = [plateauSlopes; R1.localSlope(idx1)];
    end
    if isstruct(R2) && R2.success
        idx2 = max(1, R2.i-3) : min(length(R2.localSlope), R2.j+3);
        plateauSlopes = [plateauSlopes; R2.localSlope(idx2)];
    end

    if isstruct(R1) && R1.success
        r_fit1 = R1.r(R1.i:R1.j);

        fill([r_fit1(1), r_fit1(end), r_fit1(end), r_fit1(1)], ...
             [-10, -10, 10, 10], c1, ...
             'FaceAlpha', 0.10, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        h1 = plot(R1.r, R1.localSlope, 'o-', 'Color', c1, ...
            'LineWidth', 1.2, 'MarkerSize', 3.5, 'MarkerFaceColor', c1);
        plot([r_fit1(1), r_fit1(end)], [R1.D2, R1.D2], '--', ...
            'Color', c1, 'LineWidth', 1.8, 'HandleVisibility', 'off');

        hLeg = [hLeg, h1];
        tLeg{end+1} = sprintf('%s $\\approx %.3f$', lbl1, R1.D2);
    end

    if isstruct(R2) && R2.success
        r_fit2 = R2.r(R2.i:R2.j);

        fill([r_fit2(1), r_fit2(end), r_fit2(end), r_fit2(1)], ...
             [-10, -10, 10, 10], c2, ...
             'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        h2 = plot(R2.r, R2.localSlope, 's-', 'Color', c2, ...
            'LineWidth', 1.2, 'MarkerSize', 3.5, 'MarkerFaceColor', c2);
        plot([r_fit2(1), r_fit2(end)], [R2.D2, R2.D2], '--', ...
            'Color', c2, 'LineWidth', 1.8, 'HandleVisibility', 'off');

        hLeg = [hLeg, h2];
        tLeg{end+1} = sprintf('%s $\\approx %.3f$', lbl2, R2.D2);
    end

    if ~isempty(plateauSlopes)
        minY = min(plateauSlopes);
        maxY = max(plateauSlopes);
        padding = max(0.15 * (maxY - minY), 0.05);
        ylim([max(0, minY - padding), maxY + padding]);
    end

    set(gca, 'XScale', 'log', 'FontName', 'Helvetica', 'FontSize', 8.5, 'LineWidth', 0.8);
    xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 9.5);
    ylabel('$d\log C / d\log r$', 'Interpreter', 'latex', 'FontSize', 9.5);

    grid on; box on;
    set(gca, 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.8);

    if ~isempty(hLeg)
        legend(hLeg, tLeg, 'Location', 'northwest', 'Interpreter', 'latex', ...
            'FontSize', 7.5, 'Box', 'off');
    end

    text(-0.10, 1.03, sprintf('\\textbf{(%s)}', pTag), 'Units', 'normalized', ...
        'FontSize', 10, 'Interpreter', 'latex');
end

end
