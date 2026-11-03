function ylims = drawLocalSlopePanel(R1, R2, lbl1, lbl2, ylims)
% drawLocalSlopePanel  Draw one local-slope panel (two D2 profiles, selected scaling windows and a GeoC box) in the current axes.
%
% CALLS     None
% CALLED BY plotLocalSlopeDiagnostics3Cases, makeCOVIDRegionalMainFigure
%
% INPUT
%   R1, R2     estimateD2 outputs of the lower-dimensional (blue) and higher-dimensional (red) state
%   lbl1, lbl2 TeX labels of the two D2 terms, for example 'D_2(\mu_n)'
%   ylims      [low high] y limits, or [] to set them from the two profiles
% OUTPUT
%   ylims      the y limits used
c1 = [0.12, 0.31, 0.89];
c2 = [0.85, 0.15, 0.20];
cFail = [0.50, 0.50, 0.50];
Rs = {R1, R2};
cols = {c1, c2};
marks = {'o-', 's-'};
alphas = [0.10, 0.08];
lbls = {lbl1, lbl2};

if isempty(ylims)
    s = [];
    for k = 1:2
        if hasDiagnostic(Rs{k}), v = Rs{k}.localSlope(:); s = [s; v(isfinite(v))]; end
    end
    if isempty(s)
        ylims = [0 1];
    else
        pad = max(0.15*(max(s) - min(s)), 0.05);
        ylims = [max(0, min(s) - pad), max(s) + pad];
    end
end
yLower = ylims(1); yUpper = ylims(2);

hold on;
boxLine = {'', ''};
for k = 1:2
    R = Rs{k}; col = cols{k};
    if ~hasDiagnostic(R), continue; end
    if R.success && isfinite(R.D2) && isfinite(R.i) && isfinite(R.j)
        rFit = R.r(R.i:R.j);
        fill([rFit(1), rFit(end), rFit(end), rFit(1)], [yLower, yLower, yUpper, yUpper], ...
            col, 'FaceAlpha', alphas(k), 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(R.r, R.localSlope, marks{k}, 'Color', col, 'LineWidth', 1.2, 'MarkerSize', 3.5, 'MarkerFaceColor', col);
        plot([rFit(1), rFit(end)], [R.D2, R.D2], '--', 'Color', col, 'LineWidth', 1.8);
        boxLine{k} = sprintf('{\\color[rgb]{%.3f,%.3f,%.3f} %s \\approx %.3f}', col(1), col(2), col(3), lbls{k}, R.D2);
    else
        rValid = R.r(isfinite(R.r) & R.r > 0);
        if ~isempty(rValid)
            fill([rValid(1), rValid(end), rValid(end), rValid(1)], [yLower, yLower, yUpper, yUpper], ...
                cFail, 'FaceAlpha', 0.08, 'EdgeColor', 'none');
        end
        plot(R.r, R.localSlope, marks{k}, 'Color', col, 'LineWidth', 1.2, 'MarkerSize', 3.5, 'MarkerFaceColor', col);
        finiteSlope = R.localSlope(isfinite(R.localSlope));
        if ~isempty(finiteSlope)
            yline(median(finiteSlope), ':', 'Color', col, 'LineWidth', 1.5);
        end
        boxLine{k} = sprintf('{\\color[rgb]{%.3f,%.3f,%.3f} %s = NR}', col(1), col(2), col(3), lbls{k});
    end
end

diffLine = '';
if hasDiagnostic(R1) && hasDiagnostic(R2)
    if R1.success && isfinite(R1.D2) && R2.success && isfinite(R2.D2)
        diffLine = sprintf('{\\bf\\color{black} \\Delta \\approx %.3f}', R2.D2 - R1.D2);
    else
        diffLine = '{\\bf\\color{black} \\Delta = NR}';
    end
end

boxLines = {};
if ~isempty(boxLine{2}), boxLines{end+1} = boxLine{2}; end
if ~isempty(boxLine{1}), boxLines{end+1} = boxLine{1}; end
if ~isempty(diffLine), boxLines{end+1} = diffLine; end
if ~isempty(boxLines)
    [xPos, yPos, hAlign, vAlign] = findBestTextPosition(R1, R2, yLower, yUpper);
    text(xPos, yPos, boxLines, 'Units', 'normalized', 'Interpreter', 'tex', 'FontSize', 6.5, ...
        'HorizontalAlignment', hAlign, 'VerticalAlignment', vAlign, ...
        'BackgroundColor', [1 1 1 0.25], 'EdgeColor', [0.65 0.65 0.65], 'LineWidth', 0.5, 'Margin', 1.8);
end

ylim(ylims);
ax = gca;
set(ax, 'XScale', 'log', 'FontName', 'Helvetica', 'FontSize', 8, 'LineWidth', 0.8, 'TickLabelInterpreter', 'latex');
grid on; box on;
set(ax, 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.8);
end

function tf = hasDiagnostic(R)
tf = isstruct(R) && isfield(R, 'r') && isfield(R, 'localSlope') && ~isempty(R.r) && ~isempty(R.localSlope);
end

function [xPos, yPos, hAlign, vAlign] = findBestTextPosition(R1, R2, yLower, yUpper)
rVals = []; sVals = [];
if hasDiagnostic(R1), rVals = [rVals; R1.r(:)]; sVals = [sVals; R1.localSlope(:)]; end
if hasDiagnostic(R2), rVals = [rVals; R2.r(:)]; sVals = [sVals; R2.localSlope(:)]; end
ok = isfinite(rVals) & rVals > 0 & isfinite(sVals);
rVals = rVals(ok); sVals = sVals(ok);
if isempty(rVals)
    xPos = 0.03; yPos = 0.95; hAlign = 'left'; vAlign = 'top';
    return;
end
logR = log10(rVals);
minX = min(logR); maxX = max(logR);
if maxX == minX, maxX = minX + 1; end
xNorm = (logR - minX) / (maxX - minX);
yNorm = (sVals - yLower) / max(yUpper - yLower, 1e-6);
candidates = {
    0.03, 0.95, 'left',  'top',    [0.0, 0.40], [0.60, 1.0];
    0.97, 0.95, 'right', 'top',    [0.60, 1.0], [0.60, 1.0];
    0.03, 0.05, 'left',  'bottom', [0.0, 0.40], [0.0, 0.40];
    0.97, 0.05, 'right', 'bottom', [0.60, 1.0], [0.0, 0.40]};
scores = zeros(4, 1);
for i = 1:4
    xR = candidates{i, 5}; yR = candidates{i, 6};
    scores(i) = sum(xNorm >= xR(1) & xNorm <= xR(2) & yNorm >= yR(1) & yNorm <= yR(2));
end
[~, bestIdx] = min(scores);
xPos = candidates{bestIdx, 1}; yPos = candidates{bestIdx, 2};
hAlign = candidates{bestIdx, 3}; vAlign = candidates{bestIdx, 4};
end
