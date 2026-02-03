clear all; close all; clc;
load('USACovid19s01092020f05222023.mat')
CFR_USA=Mu;    
Cases_USA_norm=(Beta - min(Beta)) ./ (max(Beta) - min(Beta));
Day_USA=Day;
clear Mu Beta Day
load('EuropeCovid19s01092020f08242025.mat')
CFR_EU=Mu;    
Cases_EU_norm=(Beta - min(Beta)) ./ (max(Beta) - min(Beta));
dates=Day;
clear Mu Beta Day
load('AsiaCovid19s01092020f08242025.mat')
 CFR_ASIA=Mu; 
 Cases_ASIA_norm=(Beta - min(Beta)) ./ (max(Beta) - min(Beta));
 Day_ASIA=Day;
 clear Mu Beta Day
load('Covid19s01092020f08242025.mat')
 CFR_Global=Mu; 
 Cases_Global_norm=(Beta - min(Beta)) ./ (max(Beta) - min(Beta));

%%


% -- Colors --
cCFR  = [0 0 0];                
cCase = [0 0.4470 0.7410];      


figW = 7.0;  
figH = 5.0;   
fs_axes = 9;  
fs_title = 10;
lnw_main = 0.5;
lnw_small = 0.5;
fontname = 'Times';  
hFig = figure('Color','w','Units','inches','Position',[1 1 figW figH],...
              'PaperPositionMode','auto','Renderer','painters');

% ==================
% TOP PANEL: Global 
% ==================
ax1 = subplot(2,3,[1 2 3],'Parent',hFig);
yyaxis(ax1,'left');
p1 = plot(ax1, dates,CFR_Global, 'LineWidth', lnw_main, 'Color', cCFR);
ylabel(ax1,'CFR (%)','FontName',fontname,'FontSize',fs_axes);

yyaxis(ax1,'right');
p2 = plot(ax1, dates, Cases_Global_norm, 'LineWidth', lnw_main, 'Color', cCase);
ylabel(ax1,'cases (scaled)','FontName',fontname,'FontSize',fs_axes);

title(ax1,'Global: Virulence and Transmission','FontWeight','bold','FontName',fontname,'FontSize',fs_title);
grid(ax1,'on');
ax1.XAxis.TickLabelFormat = 'MMM yyyy';
xlim(ax1,[min(dates) max(dates)]);

% Set axis colors explicitly
yyaxis(ax1,'left');  ax1.YAxis(1).Color = cCFR;  ax1.YLabel.Color = cCFR;
yyaxis(ax1,'right'); ax1.YAxis(2).Color = cCase; ax1.YLabel.Color = cCase;

% =====================
% BOTTOM PANEL 
% ===================
ax2 = subplot(2,3,4,'Parent',hFig);
yyaxis(ax2,'left');
plot(ax2, Day_USA, CFR_USA, 'LineWidth', lnw_small, 'Color', cCFR);
ylabel(ax2,'CFR (%)','FontName',fontname,'FontSize',fs_axes);

yyaxis(ax2,'right');
plot(ax2, Day_USA, Cases_USA_norm, 'LineWidth', lnw_small, 'Color', cCase);
ylabel(ax2,'Cases (scaled)','FontName',fontname,'FontSize',fs_axes);

title(ax2,'USA','FontName',fontname,'FontSize',fs_title);
grid(ax2,'on'); ax2.XAxis.TickLabelFormat = 'MMM yyyy';

yyaxis(ax2,'left');  ax2.YAxis(1).Color = cCFR;  ax2.YLabel.Color = cCFR;
yyaxis(ax2,'right'); ax2.YAxis(2).Color = cCase; ax2.YLabel.Color = cCase;


ax3 = subplot(2,3,5,'Parent',hFig);
yyaxis(ax3,'left');
plot(ax3, dates, CFR_EU, 'LineWidth', lnw_small, 'Color', cCFR);
ylabel(ax3,'CFR (%)','FontName',fontname,'FontSize',fs_axes);

yyaxis(ax3,'right');
plot(ax3, dates, Cases_EU_norm, 'LineWidth', lnw_small, 'Color', cCase);
ylabel(ax3,'Cases (scaled)','FontName',fontname,'FontSize',fs_axes);

title(ax3,'Europe','FontName',fontname,'FontSize',fs_title);
grid(ax3,'on'); ax3.XAxis.TickLabelFormat = 'MMM yyyy';

yyaxis(ax3,'left');  ax3.YAxis(1).Color = cCFR;  ax3.YLabel.Color = cCFR;
yyaxis(ax3,'right'); ax3.YAxis(2).Color = cCase; ax3.YLabel.Color = cCase;


ax4 = subplot(2,3,6,'Parent',hFig);
yyaxis(ax4,'left');
plot(ax4, Day_ASIA, CFR_ASIA, 'LineWidth', lnw_small, 'Color', cCFR);
ylabel(ax4,'CFR (%)','FontName',fontname,'FontSize',fs_axes);

yyaxis(ax4,'right');
plot(ax4, Day_ASIA, Cases_ASIA_norm, 'LineWidth', lnw_small, 'Color', cCase);
ylabel(ax4,'Cases (scaled)','FontName',fontname,'FontSize',fs_axes);

title(ax4,'Asia','FontName',fontname,'FontSize',fs_title);
grid(ax4,'on'); ax4.XAxis.TickLabelFormat = 'MMM yyyy';

yyaxis(ax4,'left');  ax4.YAxis(1).Color = cCFR;  ax4.YLabel.Color = cCFR;
yyaxis(ax4,'right'); ax4.YAxis(2).Color = cCase; ax4.YLabel.Color = cCase;


allAxes = [ax1 ax2 ax3 ax4];
for ax = allAxes
    set(ax, 'FontName', fontname, 'FontSize', fs_axes);
    ax.Box = 'off';
    ax.TickLabelInterpreter = 'none'; 
end

set([ax2 ax3 ax4], 'XTickLabelRotation', 20);


linkaxes(allAxes,'x');


hold(ax1,'on');
leg_h1 = plot(ax1, nan, nan, '-', 'Color', cCFR, 'LineWidth', lnw_main);
leg_h2 = plot(ax1, nan, nan, '-', 'Color', cCase, 'LineWidth', lnw_main);
legend(ax1, [leg_h1 leg_h2], {'CFR','Cases (scaled)'}, ...
       'Location','northeastoutside','FontName',fontname,'FontSize',fs_axes);
hold(ax1,'off');

drawnow;

%====================
% EXPORT 
% =====================
outbase = 'COVID_time_series';

% Vector PDF (preferred for journals)
set(hFig,'Renderer','painters'); % ensures vector output
% exportgraphics(hFig, [outbase '.pdf'], 'ContentType','vector');  % PDF (vector)
% 
% % High-res PNG  600 dpi
% exportgraphics(hFig, [outbase '.png'], 'Resolution',600);       % 600 dpi PNG
% 
% % EPS 
% print(hFig, [outbase '.eps'], '-depsc2', '-r300');
% 
% fprintf('Saved files: %s.pdf  %s.png  %s.eps\n', outbase, outbase, outbase);
