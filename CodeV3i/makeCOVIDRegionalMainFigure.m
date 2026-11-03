function fig = makeCOVIDRegionalMainFigure(results, opts)
% makeCOVIDRegionalMainFigure  Regional synthesis figure (post-mass-vaccination GeoC(mu_{n+1}|mu_n) and marginals for Asia, Europe and USA).
%
% CALLS     drawLocalSlopePanel
% CALLED BY runCOVID
%
% INPUT
%   results   struct array of analyzeCOVIDRegionLETR outputs that includes the ids 'Asia', 'Europe' and 'USA'
%   opts      struct with OutputRoot, Figures.Resolution (dpi) and Figures.ExportPDF (logical)
% OUTPUT
%   fig       figure handle (also saved to OutputRoot/MainFigures)
regionIds = {'Asia','Europe','USA'};
ids = {results.id};
fontName = 'Helvetica'; tickSize = 8; labelSize = 9;
muColor = [0.20,0.35,0.55]; betaColor = [0.70,0.35,0.20];

geoR1 = cell(1,3); geoR2 = cell(1,3); inv = cell(1,3);
allSlopes = [];
muAxisMax = 0;
for k = 1:3
    post = results(strcmp(ids, regionIds{k})).periods(3);
    geoR1{k} = post.GeoC.details.D2_mun;
    geoR2{k} = post.GeoC.details.D2_mun_munp1;
    inv{k} = post.invariantDensity;
    for R = {geoR1{k}, geoR2{k}}
        s = R{1}.localSlope(:);
        allSlopes = [allSlopes; s(isfinite(s))]; 
    end
    muAxisMax = max(muAxisMax, max(inv{k}.muVec));
end
pad = max(0.15*(max(allSlopes) - min(allSlopes)), 0.05);
ylims = [max(0, min(allSlopes) - pad), max(allSlopes) + pad];

fig = figure('Color','w','Units','centimeters','Position',[2 2 24 20],'Renderer','painters');
t = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

for c = 1:3
    nexttile;
    drawLocalSlopePanel(geoR1{c}, geoR2{c}, 'D_2(\mu_n)', 'D_2([\mu_n,\mu_{n+1}])', ylims);
    xlabel('$r$','Interpreter','latex','FontSize',9);
    if c == 1, ylabel('$d\log C / d\log r$','Interpreter','latex','FontSize',9); end
    title(regionIds{c},'Interpreter','none','FontName',fontName,'FontSize',9.5,'FontWeight','normal');
    text(-0.12,1.04,sprintf('\\textbf{(%s)}',char('a'+c-1)),'Units','normalized','FontSize',9,'Interpreter','latex');
    hold off;
end

for k = 1:3
    axM = nexttile(t,3+k);
    plot(axM,inv{k}.muVec,inv{k}.f_mu_marginal,'-','Color',muColor,'LineWidth',1.5);
    xlim(axM,[0 muAxisMax]);
    xlabel(axM,'$\mu$','Interpreter','latex','FontName',fontName,'FontSize',labelSize);
    if k == 1, ylabel(axM,'$p^*_{\mu}(\mu)$','Interpreter','latex','FontName',fontName,'FontSize',labelSize); end
    styleDensityAxes(axM,fontName,tickSize);
    text(axM,-0.12,1.04,sprintf('\\textbf{(%s)}',char('a'+3+k-1)),'Units','normalized','FontSize',9,'Interpreter','latex');
end

for k = 1:3
    axB = nexttile(t,6+k);
    plot(axB,inv{k}.betaVec,inv{k}.f_beta_marginal,'--','Color',betaColor,'LineWidth',1.5);
    xlim(axB,[0 1]);
    xlabel(axB,'$\beta$','Interpreter','latex','FontName',fontName,'FontSize',labelSize);
    if k == 1, ylabel(axB,'$p^*_{\beta}(\beta)$','Interpreter','latex','FontName',fontName,'FontSize',labelSize); end
    styleDensityAxes(axB,fontName,tickSize);
    text(axB,-0.12,1.04,sprintf('\\textbf{(%s)}',char('a'+6+k-1)),'Units','normalized','FontSize',9,'Interpreter','latex');
end

set(fig,'InvertHardcopy','off');
outDir = fullfile(opts.OutputRoot,'MainFigures');
if ~exist(outDir,'dir'), mkdir(outDir); end
exportgraphics(fig,fullfile(outDir,'FigureRegionalCOVIDMain.png'),'Resolution',opts.Figures.Resolution);
if opts.Figures.ExportPDF
    exportgraphics(fig,fullfile(outDir,'FigureRegionalCOVIDMain.pdf'),'ContentType','vector');
end
end

function styleDensityAxes(ax,fontName,tickSize)
set(ax,'FontName',fontName,'FontSize',tickSize,'LineWidth',0.7,'TickDir','out', ...
    'TickLength',[0.015 0.015],'Box','off','Layer','top','XMinorTick','on','YMinorTick','on', ...
    'TickLabelInterpreter','latex');
grid(ax,'off');
end
