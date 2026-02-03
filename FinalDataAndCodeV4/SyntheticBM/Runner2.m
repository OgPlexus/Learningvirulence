clear all; close all; clc
% Generate random data
n =10^4;
b = 1.5;
xn = rand(1,n);
yn = b*rand(1,n);
a=3.999;
%%
xnp1 = 3.99*xn.*(1-xn);
%ynp1=2*xn.*(1-xn);
ynp1=sech(xn).^2;
%%
% Iterate function
for i = 1:500
    if i ~= 1
        xn = xnp1;
        yn=ynp1;
    end
    xnp1 = 3.999*xn.*(1-xn);
    %ynp1=b*sech(xn+yn).^2;
    ynp1=b*sech(xn+yn).^2;
end
%%
%run Runner2.m or 
clear all; close all; clc;
load("XYdata1.mat")
%for plots in the article
%%
 close all
scatter(xnp1,ynp1,'.');
xlabel('$\mu$');
ylabel('$\beta$');
%title("$\mu$ vs $\beta$")
hfig = gcf;
fname = 'tradOff';
%%
close all;
scatter(xn,xnp1,50,'.');
xlabel('Virulence: $\mu_{n}$');
ylabel('Virulence: $\mu_{n+1}$');
%title("Virulance Relation")
hfig = gcf;
fname = 'muTseries';
%%
close all;
scatter(yn,ynp1,'.');
xlabel('Transmission: $\beta_{n}$');
ylabel('Transmission: $\beta_{n+1}$');
%title("Transmission relation ")
hfig = gcf;
fname = 'betaTseries';
%%
close all;
scatter3(xn,yn, ynp1,'.');
xlabel('Virulence: $\mu_{n}$');
ylb=ylabel('Transmission: $\beta_{n}$','Rotation',-30);
ylb.Position=[-0.2 0.5 0];
zlabel('Transmission: $\beta_{n+1}$');
%title("causal relationship")
hfig = gcf;
fname = 'muBetaBetap';