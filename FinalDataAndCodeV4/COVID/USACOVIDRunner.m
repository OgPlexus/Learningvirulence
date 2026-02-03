%Beta=Beta(242:end);
clear all; close all; clc;
load('USACovid19s01092020f05222023.mat')
% colors.usa    = [0.18 0.32 0.47];   
% colors.eu     = [0.20 0.60 0.60];   
% colors.asia   = [0.55 0.72 0.82];   


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
MuUsa =Mu;
BetasUSA= Betas;
%%
close all
subplot(2,1,1)
plot(Mu,'k.');xlabel('time (days)');ylabel('CFR')
subplot(2,1,2)
plot(Betas,'b.');xlabel('time (days)');ylabel('Rescaled daily new cases')
%sgtitle('Time series data');
hfig = gcf;
fname = 'USACOVIDtseriesdata1';
%% Fitting Beta=a/mu full world
% modelFun = @(a,mu)  a(3)+a(1) ./ (mu+a(2));
% 
% a0 = [2 0.1 1];
% 
% % Fit using nonlinear least squares
% a_est = nlinfit(X1, Y1s, modelFun, a0);
%%
close all
figure
scatter(X1, Y1s,20, 'filled','MarkerFaceColor',colors.usa); hold on
%mu_fit = linspace(min(X1), max(X1), 200)';
 %plot(mu_fit, modelFun(a_est, mu_fit), 'r-', 'LineWidth', 2)
xlabel('Virulence: $\mu$')
ylabel('Transmission: $\beta$')
%legend('Data',sprintf('Fitted curve: $\\beta = %.5f / \\mu$', a_est))
grid on
hfig = gcf;
fname = 'USACOVIDMuBeta';
%% Mu(n+1)=a1+a2*Mu(n)*(1-Mu(n))+a3*Beta(n)
rng(123);
X = [ X1.*(1-X1), Y1s];
y = X2; %Mu(n+1)

% LASSO with cross-validation
% [B,FitInfo] = lasso(X, y, 'CV', 5);%, 'Standardize', true);
% 
% idx = FitInfo.IndexMinMSE;
% coef = B(:, idx);               
% intercept = FitInfo.Intercept(idx);
% 
% yhat = X * coef + intercept;
%sci2latex = @(x) sprintf('%.1f\\times10^{%d}', x*10^abs(floor(log10(abs(x)))+1), floor(log10(abs(x)))+1);
% label2 = sprintf('Fitted curve: $\\mu_{n+1} \\approx 0.000015+ 1.023 \\mu_n (1-\\mu_n) - 0.0005 \\beta_n$');
% label2 = 'Fitted curve';
close all;
scatter(X1,X2,20,'filled','MarkerFaceColor',colors.usa);
xlabel('Virulence: $\mu_n$');
ylabel('Virulence: $\mu_{n+1}$');
hold on;
% plot(X1, yhat,'-','LineWidth',2);
% legend('Data',label2)

hfig = gcf;
fname = 'USACOVIDmuNmuNp1';
%%
%% Mu(n+1)=a1+a2*Mu(n)*(1-Mu(n))+a3*Beta(n)
rng(123);
X = [Y1s X1];
y = Y2s; %Mu(n+1)

% LASSO with cross-validation
% [B,FitInfo] = lasso(X, y, 'CV', 5);%, 'Standardize', true);
% 
% idx = FitInfo.IndexMinMSE;
% coef = B(:, idx);               
% intercept = FitInfo.Intercept(idx);
% 
% yhat = X * coef + intercept;
% %sci2latex = @(x) sprintf('%.1f\\times10^{%d}', x*10^abs(floor(log10(abs(x)))+1), floor(log10(abs(x)))+1);
% %label2 = sprintf('Fitted curve: $\\beta_{n+1} \\approx 0.0256+ 0.6917 \\beta_n - 0.7228 \\mu_n$');
% label2 = 'Fitted curve';
close all;
scatter(Y1s,Y2s,20,'filled','MarkerFaceColor',colors.usa);
xlabel('Transmission: $\beta_n$');
ylabel('Transmission: $\beta_{n+1}$');
hold on;
% plot(Y1s, yhat,'-','LineWidth',2);
% legend('Data',label2)

hfig = gcf;
fname = 'COVIDbetaNbetaNp1';

% close all
% scatter(X1,X2,20,'filled');
% xlabel('$\mu_{n}$');
% ylabel('$\mu_{n+1}$');
% %title("Virulance Relation")
% hfig = gcf;
% fname = 'COVIDmuNmuNp1';
%%
close all;
scatter(Y1s,Y2s,20,'filled');
xlabel('$\beta_{n}$');
ylabel('$\beta_{n+1}$');
%title("Transmission relation ")
hfig = gcf;
fname = 'COVIDbetaNbetaNp1';
%%
close all;
scatter(X1,Y1s,20,'filled');
xlabel('$\beta_{n}$');
ylabel('$\beta_{n+1}$');
%title("Transmission relation ")
hfig = gcf;
fname = 'COVIDMuNbetaN';
%%
close all;
scatter3(X1,Y1s,Y2s,'filled');
xlabel('$\mu_{n}$');
ylabel('$\beta_{n}$');
zlabel('$\beta_{n+1}$')
%title("causal relationship")
hfig = gcf;
fname = 'COVIDcausalRela3dBeta';
%%
close all;
scatter3(X1,Y1s,X2,'filled');
xlabel('$\mu_{n}$');
ylabel('$\beta_{n}$');
zlabel('$\beta_{n+1}$')
%title("causal relationship")
hfig = gcf;
fname = 'COVIDcausalRela3dMu';
%% check documentation of correlationDimension() for choosing parameter values.
dim1=correlationDimension(X1,[],1,'NumPoints',100, 'MinRadius',0.0003,'MaxRadius',0.003)
dim2=correlationDimension([X1 X2],[],1,'NumPoints',100,'MinRadius',0.0009,'MaxRadius',0.003)
GeoCxnp1Gxn=dim2-dim1
%%
dim1xy=correlationDimension([X1 Y1s],[],1,'NumPoints',100,'MinRadius',0.0015,'MaxRadius',0.003)
dim2xy=correlationDimension([X1 Y1s X2],[],1,'NumPoints',100,'MinRadius',0.001,'MaxRadius',0.01)
GeoCxnynxnp1Gxnyn=dim2xy-dim1xy
GeoCytox=GeoCxnp1Gxn-GeoCxnynxnp1Gxnyn
%%
dim1y=correlationDimension(Y1s,[],1,'NumPoints',100,'MinRadius',0.001,'MaxRadius',0.01)
dim2y=correlationDimension([Y1s Y2s],[],1,'NumPoints',100,'MinRadius',0.003,'MaxRadius',0.01)
GeoCynp1Gyn=dim2y-dim1y
%%
dim1yx=correlationDimension([X1 Y1s],[],1,'NumPoints',100,'MinRadius',0.0015,'MaxRadius',0.003)
dim2yx=correlationDimension([X1 Y1s Y2s],[],1,'NumPoints',100,'NumPoints',100,'MinRadius',0.005,'MaxRadius',0.01)
GeoCxnynynp1Gxnyn=dim2yx-dim1yx
GeoCxtoy=GeoCynp1Gyn-GeoCxnynynp1Gxnyn
%%
scatter(X1,X2)
%%
scatter(Y1s,Y2s)
%%
close all
scatter(Y1s,X2)
%%
scatter3(Y1s,X1,X2)
%%
scatter3(Y1s,X1,Y2s)
%%
A=[ones(size(X1)) X1 X1.^2 Y1];
coef=pinv(A)*X2
%%

