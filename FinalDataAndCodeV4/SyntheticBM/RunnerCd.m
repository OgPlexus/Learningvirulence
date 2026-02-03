%run Runner2.m or
clear all; close all; clc;
load("XYdata1.mat")
%%
% dim1=corelDim1(xn);
% dim2=corelDim1([xn xnp1]);
% GeoCxnp1Gxn=dim2-dim1;
dim1=correlationDimension(xn,[],1,'NumPoints',100,'MinRadius',0.00011,'MaxRadius',0.00086);
dim2=correlationDimension([xn xnp1],[],1,'NumPoints',100,'MinRadius',0.00013,'MaxRadius',0.0285);
GeoCxnp1Gxn=dim2-dim1
%%
dim1xy=correlationDimension([xn yn],[],1,'NumPoints',100,'MinRadius',0.0006,'MaxRadius',0.11);
dim2xy=correlationDimension([xn yn xnp1],[],1,'NumPoints',100,'MinRadius',0.00714,'MaxRadius',0.1299);
GeoCxnynxnp1Gxnyn=dim2xy-dim1xy
GeoCytox=GeoCxnp1Gxn-GeoCxnynxnp1Gxnyn
%%
dim1y=correlationDimension(yn,[],1,'NumPoints',100);
dim2y=correlationDimension([yn ynp1],[],1,'NumPoints',100);
GeoCynp1Gyn=dim2y-dim1y;
%%

dim1yx=correlationDimension([xn yn],[],1,'NumPoints',100,'MinRadius',0.0006,'MaxRadius',0.11);
dim2yx=correlationDimension([xn yn ynp1],[],1,'NumPoints',100,'MinRadius',0.00234,'MaxRadius',0.0605)
GeoCxnynynp1Gxnyn=dim2yx-dim1yx
GeoCxtoy=GeoCynp1Gyn-GeoCxnynynp1Gxnyn
%%
dim1=corelDim([x1 x2]);
dim2=corelDim([x1 x2 x2t]);
an1=dim2-dim1;
xlabel('x1');ylabel('x2');zlabel('x2t')
title(sprintf('$d_2(x1,x2)=%.4f, d_2(x1,x2,x2t)=%.4f$',dim1,dim2))
hfig = gcf;
fname = 'ep02Datax1x2x2t';
%%
[T H H2 H3 H4] = ComputeTransferEntropynew(xn',yn',xnp1',8)
%% Additional Figs
close all;
Kbin = 250;
data=[xn yn];
xi = linspace(min(data), max(data), Kbin);
bw1 = 10^(-2);
[f0, xi, ~] = ksdensity(data, xi, 'Bandwidth', bw1);
plot(xi,f0);hold on
histogram(data,Kbin, 'Normalization', 'pdf');
