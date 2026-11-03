% runSampleSize  Sensitivity of the directional GeoC estimates to the number of transitions N, with local-slope plots for Case 1(a).
%
% CALLS     buildSensitivitySystems, estimateGeoCTI, generateCase1Series, estimateGeoC, plotLocalSlopeAcrossN
% INPUT     COVID data files through buildSensitivitySystems (folder COVID). Settings are in the opts block below:
%           opts.Benchmarks (case ids), opts.COVID.Regions and opts.COVID.Periods (use {} for no SARS-CoV-2 records),
%           opts.Seeds (random subsamples), opts.Tau (positivity threshold), opts.NArticle (N of the article),
%           opts.NBench and opts.NCovid (values of N in the table, Inf uses every transition),
%           opts.FigureN and opts.FigureSeed (values of N and subsample seed of the Case 1(a) figure).
%           d2Bench and d2Covid hold the article settings of the dimension estimator.
% OUTPUT    ResultsSampleSize/TableSampleSize.tex, SampleSize.csv, FigureCase1aSampleSize.png and .pdf,
%           and N_<id>.mat checkpoints (an interrupted run resumes)

clear; close all; clc;
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

%%
opts.OutputRoot    = 'ResultsSampleSize';
opts.Benchmarks    = {'1a','1bi','1bii','1biii','1c'};
opts.COVID.Regions = {'Global'};
opts.COVID.Periods = {'preVacc','postVacc','full'};
opts.Seeds         = 1:3;
opts.Tau           = 0.03;
opts.NArticle      = 2500;
opts.NBench        = [250 500 1000 2500 5000];
opts.NCovid        = [200 400 800 1200 Inf];
opts.FigureN       = [250 500 1000 2500 5000 Inf];
opts.FigureSeed    = 2;

if ~exist(opts.OutputRoot, 'dir'), mkdir(opts.OutputRoot); end

d2Bench = {'Cmax',0.20,'R2min',0.98,'SlopeCVmax',0.15,'MaxSlopeDrift',0.10,'MinPlateauLength',10,'MinPairs',1000};
d2Covid = {'Cmax',0.20,'R2min',0.95,'SlopeCVmax',0.20,'MaxSlopeDrift',0.15,'MinPlateauLength',10,'MinPairs',1000};

%%
systems = buildSensitivitySystems(opts);
nSys = numel(systems);
Out = cell(nSys, 1);
sig = [numel(opts.Seeds), opts.NBench(:)', opts.NCovid(:)'];
for s = 1:nSys
    sys = systems(s);
    ckpt = fullfile(opts.OutputRoot, ['N_' sys.id '.mat']);
    if exist(ckpt, 'file')
        L = load(ckpt);
        if isfield(L.O, 'sig') && isequal(L.O.sig, sig)
            Out{s} = L.O; fprintf('[%d/%d] %s: checkpoint found.\n', s, nSys, sys.id); continue;
        end
    end
    if strcmp(sys.group, 'covid'), Nv = opts.NCovid; d2 = d2Covid; else, Nv = opts.NBench; d2 = d2Bench; end
    R = nan(numel(Nv), numel(opts.Seeds), 2);
    t0 = tic;
    for v = 1:numel(Nv)
        for k = 1:numel(opts.Seeds)
            R(v, k, :) = runOne(sys.Mu, sys.Beta, Nv(v), opts.Seeds(k), d2);
        end
    end
    fprintf('[%d/%d] %s done (%.0f s)\n', s, nSys, sys.label, toc(t0));
    O = struct('id', sys.id, 'label', sys.label, 'group', sys.group, 'Nv', Nv, ...
        'nTrans', numel(sys.Mu) - 1, 'R', R, 'sig', sig);
    save(ckpt, 'O'); Out{s} = O;
end

writeSampleTable(Out, opts, fullfile(opts.OutputRoot, 'TableSampleSize.tex'));
exportCSV(Out, fullfile(opts.OutputRoot, 'SampleSize.csv'));

%%
[Mu, Beta] = generateCase1Series('1a');
nTrans = numel(Mu) - 1;
nFig = numel(opts.FigureN);
detList = cell(1, nFig);
Nused = zeros(1, nFig);
ws = warning('off', 'conditionalD2:noPlateau');
for q = 1:nFig
    [~, detList{q}] = estimateGeoC(Mu, Beta, 'Subsample', opts.FigureN(q), 'Seed', opts.FigureSeed, d2Bench{:});
    Nused(q) = min(opts.FigureN(q), nTrans);
end
warning(ws);

fig = plotLocalSlopeAcrossN(detList, Nused);
exportgraphics(fig, fullfile(opts.OutputRoot, 'FigureCase1aSampleSize.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(opts.OutputRoot, 'FigureCase1aSampleSize.pdf'), 'ContentType', 'vector');
close(fig);

%%
function v = runOne(Mu, Beta, N, seed, d2)
try
    o = estimateGeoCTI(Mu, Beta, 'Subsample', N, 'Seed', seed, d2{:});
    v = [o.B2M o.M2B];
catch
    v = [NaN NaN];
end
end

function s = texEsc(s)
s = strrep(s, '&', '\&'); s = strrep(s, '%', '\%'); s = strrep(s, '_', '\_');
end

function s = cellTx(x)
if isnan(x), s = 'n.p.'; else, s = sprintf('%.2f', x); end
end

function s = listN(Nv)
s = strjoin(arrayfun(@(v) sprintf('%d', v), Nv, 'UniformOutput', false), ', ');
end

function s = classChange(x, skip, iB, tau, Nv)
x = x(:)'; skip = skip(:)'; Nv = Nv(:)';
if isempty(iB) || isnan(x(iB)), s = 'n.a.'; return; end
pos = x > tau;
chg = ~skip & isfinite(x) & (pos ~= pos(iB));
np  = ~skip & isnan(x);
parts = {};
if any(chg)
    if pos(iB), a = '$+\to0$'; else, a = '$0\to+$'; end
    parts{end+1} = [a ' (' listN(Nv(chg)) ')'];
end
if any(np)
    parts{end+1} = ['n.p. (' listN(Nv(np)) ')'];
end
if isempty(parts), s = 'None'; else, s = strjoin(parts, ', '); end
end

function writeSampleTable(Out, opts, file)
grp = {'benchmark', 'covid'};
panelTitle = {'A. Benchmarks', 'B. SARS-CoV-2'};
L = {};
L{end+1} = '\begin{table}[htbp]';
L{end+1} = '\centering';
L{end+1} = '\footnotesize';
L{end+1} = '\setlength{\tabcolsep}{4pt}';
L{end+1} = ['\caption{Sensitivity of the directional GeoC estimates to the number of transitions $N$ used in each computation. ' ...
    'Only $N$ is varied and all other settings are those in \Cref{SItab:numerical_settings}. ' ...
    'The value used in the article is marked with an asterisk, and ``All'''' uses every transition of the record. ' ...
    'Entries are means over ' num2str(numel(opts.Seeds)) ' random subsamples. ' ...
    '``n.p.'''' means that no stable scaling region was found, and for subsampled entries that applies when it occurred in at least one subsample. ' ...
    'A dot indicates that $N$ exceeds the record length, so the entry equals the one for the full record. ' ...
    'The last column reports the change of class relative to the article setting, where $+$ denotes a value above ' num2str(opts.Tau) ...
    ' and $0$ a value of ' num2str(opts.Tau) ' or below, including negative values. ' ...
    'The values of $N$ at which the change occurs are given in parentheses. ' ...
    '``None'''' means that the class is the same at every value of $N$, and an entry of n.p. means that the class is unchanged at every value with a plateau.}'];
L{end+1} = '\label{SItab:sampleSize}';
first = true;
for g = 1:2
    idx = find(cellfun(@(o) strcmp(o.group, grp{g}), Out));
    if isempty(idx), continue; end
    if g == 1, Nv = opts.NBench; iB = find(Nv == opts.NArticle, 1); else, Nv = opts.NCovid; iB = find(isinf(Nv), 1); end
    nV = numel(Nv);
    hdr = arrayfun(@(x) sprintf('%d', x), Nv, 'UniformOutput', false);
    hdr(isinf(Nv)) = {'All'};
    if ~isempty(iB), hdr{iB} = [hdr{iB} '$^{*}$']; end
    if ~first, L{end+1} = '\vspace{6pt}'; end
    first = false;
    L{end+1} = ['\begin{tabular}{@{}ll' repmat('c', 1, nV) '>{\raggedright\arraybackslash}p{2.2cm}@{}}'];
    L{end+1} = '\toprule';
    L{end+1} = ['\multicolumn{' num2str(nV + 3) '}{@{}l}{\textbf{' panelTitle{g} ', number of transitions $N$}} \\'];
    L{end+1} = '\midrule';
    L{end+1} = ['System & Flow & ' strjoin(hdr, ' & ') ' & Class change \\'];
    L{end+1} = '\midrule';
    for q = 1:numel(idx)
        O = Out{idx(q)};
        capped = isfinite(Nv) & Nv >= O.nTrans;
        if ~isempty(iB), capped(iB) = false; end
        anyNaN = squeeze(any(isnan(O.R), 2));
        M = squeeze(mean(O.R, 2, 'omitnan')); M(anyNaN) = NaN;
        for d = 1:2
            if d == 1, nm = texEsc(O.label); fl = '$\beta\to\mu$'; else, nm = ''; fl = '$\mu\to\beta$'; end
            cells = arrayfun(@cellTx, M(:, d), 'UniformOutput', false);
            cells(capped) = {'$\cdot$'};
            cc = classChange(M(:, d), capped, iB, opts.Tau, Nv);
            L{end+1} = [nm ' & ' fl ' & ' strjoin(cells(:)', ' & ') ' & ' cc ' \\'];
        end
        if q < numel(idx), L{end+1} = '\addlinespace[1pt]'; end
    end
    L{end+1} = '\bottomrule';
    L{end+1} = '\end{tabular}';
end
L{end+1} = '\end{table}';
fid = fopen(file, 'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
end

function exportCSV(Out, file)
rows = {};
for s = 1:numel(Out)
    O = Out{s};
    for v = 1:numel(O.Nv)
        for d = 1:2
            x = squeeze(O.R(v, :, d));
            rows(end+1, :) = {O.id, O.Nv(v), d, mean(x, 'omitnan'), std(x, 'omitnan'), nnz(isnan(x))};
        end
    end
end
T = cell2table(rows, 'VariableNames', {'System','N','Direction_1b2m_2m2b','Mean','SD','nNoPlateau'});
writetable(T, file);
end
