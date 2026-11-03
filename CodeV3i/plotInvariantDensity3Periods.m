function fig = plotInvariantDensity3Periods(result)
% plotInvariantDensity3Periods  3-by-3 invariant-density figure of one region (rows joint, virulence and transmission marginals; columns full, pre- and post-mass-vaccination).
%
% CALLS     batlow, addPanelLabel
% CALLED BY makeCOVIDRegionFigures
%
% INPUT
%   result   struct from analyzeCOVIDRegionLETR with raw.Mu and periods(1:3).label, periods(1:3).invariantDensity
%            (estimateInvariantDensity2D output)
% OUTPUT
%   fig      figure handle (requires MATLAB R2022a or later for clim)
periods = result.periods;
fig = figure('Color','w','Units','centimeters','Position',[1.5 1.5 18 18],'Renderer','painters');
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
fontName = 'Arial'; tickSize = 8.0; labelSize = 9.5;
muColor = [0.20 0.35 0.55]; betaColor = [0.70 0.35 0.20];
muMax = max(result.raw.Mu,[],'omitnan');
betaMax = 1;
allDensity = [];
for p = 1:3
    Z = periods(p).invariantDensity.invariantDensity;
    allDensity = [allDensity; Z(isfinite(Z))];
end
jointColorMax = max(allDensity);
colormap(fig,batlow(256));

for p = 1:3
    invOut = periods(p).invariantDensity;

    ax1 = nexttile(tl,p);
    contourf(ax1,invOut.MU,invOut.BETA,invOut.invariantDensity,25,'LineColor','none');
    xlim(ax1,[0 muMax]); ylim(ax1,[0 betaMax]); clim(ax1,[0 jointColorMax]);
    xlabel(ax1,'$\mu$','Interpreter','latex','FontName',fontName,'FontSize',labelSize);
    ylabel(ax1,'$\beta$','Interpreter','latex','FontName',fontName,'FontSize',labelSize);
    title(ax1,periods(p).label,'Interpreter','none','FontName',fontName,'FontSize',9.5,'FontWeight','normal');
    formatDensityAxes(ax1,fontName,tickSize);
    if p == 3
        cb = colorbar(ax1,'eastoutside');
        cb.Label.String = '$p^*(\mu,\beta)$'; cb.Label.Interpreter = 'latex';
        cb.Label.FontName = fontName; cb.Label.FontSize = 8.5;
        cb.FontName = fontName; cb.FontSize = 7.0; cb.LineWidth = 0.6; cb.TickDirection = 'out';
    end
    addPanelLabel(ax1,sprintf('(%s)',char('a' + p - 1)));

    ax2 = nexttile(tl,3+p);
    plot(ax2,invOut.muVec,invOut.f_mu_marginal,'-','Color',muColor,'LineWidth',1.65);
    xlim(ax2,[0 muMax]);
    xlabel(ax2,'$\mu$','Interpreter','latex','FontName',fontName,'FontSize',labelSize);
    ylabel(ax2,'$p^*_{\mu}(\mu)$','Interpreter','latex','FontName',fontName,'FontSize',labelSize);
    formatDensityAxes(ax2,fontName,tickSize);
    addPanelLabel(ax2,sprintf('(%s)',char('a' + 3 + p - 1)));

    ax3 = nexttile(tl,6+p);
    plot(ax3,invOut.betaVec,invOut.f_beta_marginal,'--','Color',betaColor,'LineWidth',1.65);
    xlim(ax3,[0 betaMax]);
    xlabel(ax3,'$\beta$','Interpreter','latex','FontName',fontName,'FontSize',labelSize);
    ylabel(ax3,'$p^*_{\beta}(\beta)$','Interpreter','latex','FontName',fontName,'FontSize',labelSize);
    formatDensityAxes(ax3,fontName,tickSize);
    addPanelLabel(ax3,sprintf('(%s)',char('a' + 6 + p - 1)));
end
drawnow;
end

function formatDensityAxes(ax,fontName,tickSize)
set(ax,'FontName',fontName,'FontSize',tickSize,'LineWidth',0.75,'TickDir','out', ...
    'TickLength',[0.015 0.015],'Box','off','Layer','top','XMinorTick','on','YMinorTick','on');
grid(ax,'off');
pbaspect(ax,[1 1 1]);
end
