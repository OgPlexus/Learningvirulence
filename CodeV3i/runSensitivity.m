% runSensitivity  One-at-a-time sensitivity of the directional GeoC estimates to Cmax and to the local-slope half-window.
%
% CALLS     buildSensitivitySystems, estimateGeoCTI
% INPUT     COVID data files through buildSensitivitySystems (folder COVID). Settings are in the opts block below:
%           opts.Benchmarks (case ids), opts.COVID.Regions and opts.COVID.Periods (use {} for no SARS-CoV-2 records),
%           opts.Seeds (random subsamples, benchmarks only), opts.Tau (positivity threshold), and the sweep list Sw.
%           baseBench and baseCovid hold the article settings [Sub Cmax R2min SlopeCV Drift Theiler NumRadii HalfWin RadiusK].
% OUTPUT    ResultsOneAtATime/TableSweep_Cmax.tex, TableSweep_HalfWindow.tex (LaTeX tables), FigOneAtATime.pdf and .png,
%           OneAtATime.csv, and OAT_<id>.mat checkpoints (an interrupted run resumes)

clear; close all; clc;

%%
opts.OutputRoot    = 'ResultsOneAtATime';
opts.Benchmarks    = {'1a','1bi','1bii','1biii','1c'};
opts.COVID.Regions = {'Global'};
opts.COVID.Periods = {'preVacc','postVacc','full'};
opts.Seeds         = 1:3;
opts.Tau           = 0.03;

if ~exist(opts.OutputRoot,'dir'), mkdir(opts.OutputRoot); end

baseBench = [2500 0.20 0.98 0.15 0.10 0 80 2 10];
baseCovid = [2500 0.20 0.95 0.20 0.15 0 80 2 10];

Sw = struct('name',{},'col',{},'vals',{},'xlab',{},'tex',{});
Sw(1) = struct('name','Cmax','col',2,'vals',[0.05 0.10 0.15 0.20 0.30 0.40], ...
    'xlab','C_{max}','tex','$C_{\max}$');
Sw(2) = struct('name','HalfWindow','col',8,'vals',[1 2 3 4 5], ...
    'xlab','Local-slope half-window','tex','local-slope half-window');
nSw = numel(Sw);

%%
systems = buildSensitivitySystems(opts);
nSys = numel(systems);
Out = cell(nSys,1);
sig = [nSw, numel(opts.Seeds), cellfun(@numel,{Sw.vals}), [Sw.col]];
for s = 1:nSys
    sys = systems(s);
    ckpt = fullfile(opts.OutputRoot, ['OAT_' sys.id '.mat']);
    if exist(ckpt,'file')
        L = load(ckpt);
        if isfield(L.O,'sig') && isequal(L.O.sig, sig)
            Out{s} = L.O; fprintf('[%d/%d] %s: checkpoint found.\n', s, nSys, sys.id); continue;
        end
    end
    if strcmp(sys.group,'covid'), base = baseCovid; seeds = 1; else, base = baseBench; seeds = opts.Seeds; end
    R = cell(nSw,1);
    t0 = tic;
    for w = 1:nSw
        nV = numel(Sw(w).vals);
        R{w} = nan(nV, numel(seeds), 2);
        for v = 1:nV
            st = base; st(Sw(w).col) = Sw(w).vals(v);
            for k = 1:numel(seeds)
                R{w}(v,k,:) = runOne(sys.Mu, sys.Beta, st, seeds(k));
            end
        end
    end
    fprintf('[%d/%d] %s done (%.0f s)\n', s, nSys, sys.label, toc(t0));
    O = struct('id',sys.id,'label',sys.label,'group',sys.group,'R',{R},'sig',sig);
    save(ckpt,'O'); Out{s} = O;
end

%%
for w = 1:nSw
    writeSweepTable(Out, Sw(w), baseBench(Sw(w).col), w, opts, ...
        fullfile(opts.OutputRoot, ['TableSweep_' Sw(w).name '.tex']));
end
exportCSV(Out, Sw, fullfile(opts.OutputRoot,'OneAtATime.csv'));
makeFigure(Out, Sw, baseBench, opts);

%%
function v = runOne(Mu, Beta, s, seed)
try
    o = estimateGeoCTI(Mu, Beta, 'Subsample',s(1),'Seed',seed,'Cmax',s(2), ...
        'R2min',s(3),'SlopeCVmax',s(4),'MaxSlopeDrift',s(5),'TheilerWindow',s(6), ...
        'NumRadii',s(7),'LocalSlopeHalfWindow',s(8),'RadiusK',s(9), ...
        'MinPlateauLength',10,'MinPairs',1000);
    v = [o.B2M o.M2B];
catch
    v = [NaN NaN];
end
end

function s = texEsc(s)
s = strrep(s,'&','\&'); s = strrep(s,'%','\%'); s = strrep(s,'_','\_');
end

function s = cellTx(x)
if isnan(x), s = 'n.p.'; else, s = sprintf('%.2f', x); end
end

function writeSweepTable(Out, sw, baseVal, w, opts, file)
nV = numel(sw.vals);
valTxt = arrayfun(@(x) sprintf('%g',x), sw.vals, 'UniformOutput', false);
hdr = valTxt;
iB = find(abs(sw.vals - baseVal) < 1e-12, 1);
if ~isempty(iB), hdr{iB} = [hdr{iB} '$^{*}$']; end
L = {};
L{end+1} = '\begin{table}[htbp]';
L{end+1} = '\centering';
L{end+1} = '\footnotesize';
L{end+1} = '\setlength{\tabcolsep}{4pt}';
L{end+1} = ['\caption{Sensitivity of the directional GeoC estimates to the ' sw.tex '. Only this setting is varied and all others ' ...
    'are those of the article. The value used in the article is marked with an asterisk. Entries are means over ' ...
    num2str(numel(opts.Seeds)) ' random subsamples for the benchmarks and single runs on the full series for SARS-CoV-2. ' ...
    '``n.p.'''' means that no stable scaling region was found for the setting (in at least one subsample for the benchmarks). ' ...
    'The last column is the largest absolute difference from the article value among the settings with a plateau.}'];
L{end+1} = ['\label{SItab:oat_' sw.name '}'];
L{end+1} = ['\begin{tabular}{ll' repmat('c',1,nV) 'c}'];
L{end+1} = '\toprule';
L{end+1} = ['System & Flow & ' strjoin(hdr,' & ') ' & Max change \\'];
L{end+1} = '\midrule';
for s = 1:numel(Out)
    R = Out{s}.R{w};
    anyNaN = squeeze(any(isnan(R),2));
    M = squeeze(mean(R,2,'omitnan')); M(anyNaN) = NaN;
    for d = 1:2
        if d == 1, nm = texEsc(Out{s}.label); fl = '$\beta\to\mu$'; else, nm = ''; fl = '$\mu\to\beta$'; end
        cells = arrayfun(@cellTx, M(:,d), 'UniformOutput', false);
        if ~isempty(iB) && isfinite(M(iB,d))
            mc = sprintf('%.2f', max(abs(M(isfinite(M(:,d)),d) - M(iB,d))));
        else
            mc = '--';
        end
        L{end+1} = [nm ' & ' fl ' & ' strjoin(cells',' & ') ' & ' mc ' \\']; %#ok<AGROW>
    end
    if s < numel(Out), L{end+1} = '\addlinespace[2pt]'; end
end
L{end+1} = '\bottomrule';
L{end+1} = '\end{tabular}';
L{end+1} = '\end{table}';
fid = fopen(file,'w'); fprintf(fid,'%s\n',strjoin(L,newline)); fclose(fid);
end

function exportCSV(Out, Sw, file)
rows = {};
for s = 1:numel(Out)
    for w = 1:numel(Sw)
        R = Out{s}.R{w};
        for v = 1:numel(Sw(w).vals)
            for d = 1:2
                x = squeeze(R(v,:,d));
                rows(end+1,:) = {Out{s}.id, Sw(w).name, Sw(w).vals(v), d, ...
                    mean(x,'omitnan'), std(x,'omitnan'), nnz(isnan(x))};
            end
        end
    end
end
T = cell2table(rows,'VariableNames', ...
    {'System','Parameter','Value','Direction_1b2m_2m2b','Mean','SD','nNoPlateau'});
writetable(T,file);
end

function makeFigure(Out, Sw, baseBench, opts)
nSw = numel(Sw);
grp = {'benchmark','covid'};
dirName = {'GeoC_{\beta\rightarrow\mu}','GeoC_{\mu\rightarrow\beta}'};
col = [0.00 0.45 0.70; 0.84 0.37 0.00; 0.00 0.62 0.45; 0.80 0.47 0.65; 0.35 0.35 0.35];
mk = {'o','s','^','d','v'};
fig = figure('Color','w','Units','centimeters','Position',[2 2 24 5.6*nSw]);
tl = tiledlayout(nSw,4,'TileSpacing','compact','Padding','compact');
for w = 1:nSw
    for g = 1:2
        idx = find(cellfun(@(o) strcmp(o.group,grp{g}), Out));
        for d = 1:2
            ax = nexttile; hold(ax,'on');
            for q = 1:numel(idx)
                O = Out{idx(q)}; R = O.R{w};
                M = squeeze(mean(R,2,'omitnan')); S = squeeze(std(R,0,2,'omitnan'));
                anyNaN = squeeze(any(isnan(R),2)); M(anyNaN) = NaN;
                c = col(mod(q-1,size(col,1))+1,:);
                if g == 1
                    errorbar(ax, Sw(w).vals, M(:,d), S(:,d), ['-' mk{q}], 'Color',c, ...
                        'MarkerFaceColor',c,'MarkerSize',4,'LineWidth',1,'CapSize',2, ...
                        'DisplayName',O.label);
                else
                    plot(ax, Sw(w).vals, M(:,d), ['-' mk{q}], 'Color',c, ...
                        'MarkerFaceColor',c,'MarkerSize',4,'LineWidth',1,'DisplayName',O.label);
                end
            end
            yline(ax,0,'-','Color',[0.7 0.7 0.7],'HandleVisibility','off');
            yline(ax,opts.Tau,':','Color',[0.5 0.5 0.5],'HandleVisibility','off');
            xline(ax,baseBench(Sw(w).col),'--','Color',[0.4 0.4 0.4],'HandleVisibility','off');
            set(ax,'FontSize',8,'Box','off');
            xlim(ax,[min(Sw(w).vals) max(Sw(w).vals)] + 0.04*(max(Sw(w).vals)-min(Sw(w).vals))*[-1 1]);
            if w == nSw, xlabel(ax,Sw(w).xlab,'FontSize',8); end
            if g == 1 && d == 1, ylabel(ax,'Directional GeoC','FontSize',8); end
            if w == 1
                if g == 1, gl = 'Benchmarks'; else, gl = 'SARS-CoV-2'; end
                title(ax,[gl ', ' dirName{d}],'FontSize',8,'FontWeight','normal');
            end
            if w == 1 && d == 2
                lg = legend(ax,'Location','best','FontSize',6,'Box','off');
            end
        end
    end
end
title(tl,['Dashed vertical line: article value. Dotted horizontal line: classification threshold ' ...
    num2str(opts.Tau) '. Gaps: no stable plateau.'],'FontSize',8,'FontWeight','normal');
print(fig, fullfile(opts.OutputRoot,'FigOneAtATime.pdf'), '-dpdf', '-painters');
print(fig, fullfile(opts.OutputRoot,'FigOneAtATime.png'), '-dpng', '-r300');
end
