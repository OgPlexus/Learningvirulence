% 1D data
Kbin = 200;
%x = xn;
x = yn;
xi = linspace(min(x), max(x), Kbin);
bw1 = 10^(-2);

[density, xi] = ksdensity(x,xi, 'Bandwidth', bw1);

density_at_x = ksdensity(x, x);
%%
close all

parula_cmap = parula(256); 
inverse_parula_cmap = flipud(parula_cmap); 
figure;
scatter(x, zeros(size(x))+0.1, 1000, density_at_x, 'Marker','|'); 
colormap("parula"); %clim([-0.1,2])
cb=colorbar;
cb.Label.String = 'Density';
%title('Density Colormap Plot');
%xlabel('Virulence: $\mu$');
xlabel('Transmission: $\beta$');
hold on
plot(xi,density,'LineWidth',2)
ylabel('Density');
legend('data','estimated density','Location','best')
hfig = gcf;
%fname = 'DensityMun';
fname = 'Densitybetaun';
%%

ax = gca; 
ax.YColor = 'none'; 
set(gca, 'YTick', []); 
%%

set(gca, 'Position', [0.1 0.2 0.8 0.1]); 
%% 2d data
x = xn;
y = yn;
[density, xi, yi] = ksdensity([x, y], [x, y]);
%%
close all
figure;
scatter(x, y, 10, density, 'filled'); 
colormap("parula"); %clim([-0.1,2])
cb=colorbar;
cb.Label.String = 'Density';
%title('Density Colormap Plot');
xlabel('Virulence: $\mu$');
ylabel('Transmission: $\beta$');
hfig = gcf;
fname = 'DensityMunBetan';