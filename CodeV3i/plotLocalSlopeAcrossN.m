function fig = plotLocalSlopeAcrossN(detailsList, Nvals)
% plotLocalSlopeAcrossN  Local-slope comparison of D2(beta_n) and D2([beta_n,beta_{n+1}]) for several numbers of transitions N.
%
% CALLS     None
% CALLED BY runSampleSize
%
% INPUT
%   detailsList   cell array of estimateGeoC details structs, one per value of N
%   Nvals         numbers of transitions used, one per element of detailsList
% OUTPUT
%   fig           figure handle, one panel per N with the common y limits of all panels

nN = numel(detailsList);
nc = min(3, nN);
nr = ceil(nN / nc);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 6*nc, 5.4*nr]);
tiledlayout(nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');

c1 = [0.12, 0.31, 0.89];
c2 = [0.85, 0.15, 0.20];
lbl1 = '$D_2(\beta_n)$';
lbl2 = '$D_2([\beta_n,\beta_{n+1}])$';

allSlopes = [];
for q = 1:nN
    allSlopes = [allSlopes; windowSlopes(detailsList{q}.D2_betan); windowSlopes(detailsList{q}.D2_betan_betanp1)];
end
if isempty(allSlopes)
    yl = [];
else
    padding = max(0.15 * (max(allSlopes) - min(allSlopes)), 0.05);
    yl = [max(0, min(allSlopes) - padding), max(allSlopes) + padding];
end

for q = 1:nN
    nexttile;
    hold on;

    R1 = detailsList{q}.D2_betan;
    R2 = detailsList{q}.D2_betan_betanp1;
    hLeg = [];
    tLeg = {};

    if R1.success
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

    if R2.success
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

    if ~isempty(yl), ylim(yl); end

    set(gca, 'XScale', 'log', 'FontName', 'Helvetica', 'FontSize', 8.5, 'LineWidth', 0.8);
    xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 9.5);
    if mod(q - 1, nc) == 0
        ylabel('$d\log C / d\log r$', 'Interpreter', 'latex', 'FontSize', 9.5);
    end
    grid on; box on;
    set(gca, 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.8);

    if ~isempty(hLeg)
        legend(hLeg, tLeg, 'Location', 'best', 'Interpreter', 'latex', ...
            'FontSize', 6.5, 'Box', 'off');
    end

    if R1.success && R2.success
        ttl = sprintf('$N=%d$, $\\Delta\\approx %.3f$', Nvals(q), R2.D2 - R1.D2);
    else
        ttl = sprintf('$N=%d$, $\\Delta$ not resolved', Nvals(q));
        note = {};
        if ~R1.success, note{end+1} = [lbl1 ' no plateau']; end
        if ~R2.success, note{end+1} = [lbl2 ' no plateau']; end
        text(0.97, 0.05, note, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', 'FontSize', 6.5, 'Interpreter', 'latex');
    end
    title(ttl, 'Interpreter', 'latex', 'FontSize', 9, 'FontWeight', 'normal');

    text(-0.10, 1.03, sprintf('\\textbf{(%s)}', char('a' + q - 1)), 'Units', 'normalized', ...
        'FontSize', 10, 'Interpreter', 'latex');
end

end

function s = windowSlopes(R)
s = [];
if isstruct(R) && R.success
    idx = max(1, R.i - 3) : min(length(R.localSlope), R.j + 3);
    s = R.localSlope(idx);
    s = s(:);
end
end
