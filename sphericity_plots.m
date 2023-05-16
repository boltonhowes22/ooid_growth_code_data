% Plot the sphericity of ooids from Bahamas, GSL, and Turks and Caicos
clear; close all; clc

% Andros Island
andros = readtable('data/camsizer/andros_camsizer.xlsx');
andros_spher = table2array(andros(:,5));

% Turks and Caicos
tc = readtable('data/camsizer/Trower2018_dataSupplement.xlsx'); % Trower 2018
tc_all = load('data/camsizer/tc_all.mat'); % From Jamison-Todd 2020
tc_spher = [table2array(tc(:,6)); tc_all.sample'];

% Great Salt Lake
gsl = readtable('data/camsizer/gsl_camsizer_data.csv');
gsl_spher = gsl.Var6;

%%
figure(1); clf
box on;hold on; grid on 

% Format the bins
lims = [min([tc_spher; andros_spher]),max([tc_spher; andros_spher])];
range = lims(2)-lims(1);
edges = linspace(lims(1)-0.1*range,lims(2)+0.1*range,100);
centers = movmean(edges,2,'endpoints','discard');

% Histograms
hc_andros = histcounts(andros_spher,edges);
hc_tc = histcounts(tc_spher,edges);
hc_gsl = histcounts(gsl_spher,edges);
xinterp = linspace(edges(1),edges(end),200);


k = 15;
y_andros = interp1(centers,movmean(hc_andros,k,'endpoints','shrink'),xinterp);
y_tc = interp1(centers,movmean(hc_tc,k,'endpoints','shrink'),xinterp);
y_gsl = interp1(centers,movmean(hc_gsl,k,'endpoints','shrink'),xinterp);

% Remove the NaNs
y_andros(isnan(y_andros)) = 0;
y_tc(isnan(y_tc)) = 0;
y_gsl(isnan(y_gsl)) = 0;

% The colors
c1 = hex2rgb('#0057ba');
c2 = hex2rgb('#65a98f');
c3 = hex2rgb('#a93400');

% Plot the distributions
gwin = gausswin(10); gwin = gwin/sum(gwin);
pA = fill(xinterp,conv(y_andros/nansum(y_andros),gwin,'same'),c1,'facealpha',0.2,'edgecolor',c1);
pt = fill(xinterp,conv(y_tc/nansum(y_tc),gwin,'same'),c2,'facealpha',0.2,'edgecolor',c2);
pg = fill(xinterp,conv(y_gsl/nansum(y_gsl),gwin,'same'),c3,'facealpha',0.2,'edgecolor',c3);

% Cube and Sphere
line([.806 .806], [0 1], 'color', 'k')
line([1 1], [0 1], 'color', 'k')

% Format
xlabel('Sphericity','interpreter','none');
ylim([0 .04])
xlim([.75 1.02])

set(gca, 'YTick', [], 'fontsize', 14)
pbaspect([1 1 1])