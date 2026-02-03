clear all; close all; clc;
load('Covid19s01092020f08242025.mat')

X1=Mu(1:end-1);
X2=Mu(2:end);
Y1=Beta(1:end-1);
Y2=Beta(2:end);
Y1_min = min(Y1);
Y1_max = max(Y1);
Y1s = (Y1 - Y1_min) ./ (Y1_max - Y1_min);

Y2_min = min(Y2);
Y2_max = max(Y2);
Y2s = (Y2 - Y2_min) ./ (Y2_max - Y2_min);
Beta_min=min(Beta);Beta_max=max(Beta);
Betas=(Beta - Beta_min) ./ (Beta_max - Beta_min);
%%

nbins = 100;
Va=Mu;
Va1=X1;
Va2=X2;
% Va=Betas;
% Va1=Y1s;
% Va2=Y2s;

edges   = linspace(min(Va), max(Va), nbins+1);
centers = (edges(1:end-1) + edges(2:end)) / 2; 

i1 = discretize(Va1, edges);
i2 = discretize(Va2, edges);
C  = accumarray([i1 i2], 1, [nbins nbins]);
P  = C ./ max(1, sum(C,2));   

% invariant density
[V,D] = eig(P');                       
[~,idx] = min(abs(diag(D)-1));         
pi = real(V(:,idx));
pi = pi / sum(pi);                     % normalize to sum 1
%%
close all
figure;
bar(centers, pi, 'FaceColor',[0.3 0.3 0.8], 'EdgeColor','none');
xlabel('Virulence: $\mu$','FontSize',12);
%xlabel('Transmission: $\beta$','FontSize',12);
ylabel('Invariant density','FontSize',12);
%title('Estimated Invariant Density');
hfig = gcf;
fname = 'COVIDMuInvDen';
%fname = 'COVIDBetaInvDen';