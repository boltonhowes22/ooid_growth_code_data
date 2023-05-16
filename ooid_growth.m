clear; close all;

colors = brewermap([9], 'Accent');
colors = colors([1 2 3 5 7 8],:);   

%% 1. load verts
load('data/reconstructions/gatesburg_verts.mat');
load('data/reconstructions/etina_verts.mat');
load('data/reconstructions/ethiopia_verts.mat');
load('data/reconstructions/triassic_verts.mat');
load('data/reconstructions/gsl_verts.mat');

L = 7;
ind_c = []; % Indeces for coloring the plot


all_verts = [gatesburg_verts,triassic_verts, etina_verts, ethiopia_verts, gsl_verts];
%all_verts = [etina_verts(1:4)];

% This function does all the SPHARM and Forecasting. 
[all_rot_inv , all_grown_rot_inv, all_lmcosi, all_grown_lmcosi] = vertexNormalForecastSPHARM(all_verts, L);


%% Get the indeces of the different samples

ind_col = [1 sum(cellfun(@numel,gatesburg_verts))-length(gatesburg_verts)];
ind_col(3) = ind_col(2) + sum(cellfun(@numel,triassic_verts))-length(triassic_verts);
ind_col(4) = ind_col(3) + sum(cellfun(@numel,etina_verts))-length(etina_verts);
ind_col(5) = ind_col(4) + sum(cellfun(@numel,ethiopia_verts))-length(ethiopia_verts);
ind_col(6) = ind_col(5) + sum(cellfun(@numel,gsl_verts))-length(gsl_verts);


%% Get the variables for plotting 

for order = 0:4
    X = order +1;
    [b{X}, delta_store{X}, y_previous{X}, y{X}, y_pred{X}, increment{X}] = extractPlotVariables(all_verts, all_rot_inv, all_grown_rot_inv, order);

    rmse_all(X) = RMSE(y{X},y_pred{X});
    mape_all(X) = MAPE(y{X},y_pred{X});
end

%% Make L0 Figure 
order = 0;
X = order +1;

% Threshold percentile for largest L0
pt_thresh = 75;

% isolate values that are below the L0 threshold
ind = y_previous{X}<prctile(y_previous{X}, pt_thresh);;


%% Plot Order 0
order = 0;

figure(1); clf
hold on; box on; grid on;

[mdl, ~, ~, ~] = special_plotter(b, delta_store, increment, ind_col, ind, colors, order);

% Formatting 
xlim([0 10])
ylim([-.3 .3])
xlabel([])
set(gca, 'YTick', [-.2 0 ,.2], 'fontsize', 12)

% R^2 and P-value
text(8, .2, sprintf('R^2: %0.3f',mdl.Rsquared.Ordinary), 'fontsize', 12)
text(8, .145, sprintf('p-value: %0.3f',mdl.Coefficients.pValue(2)), 'fontsize', 12)

pbaspect([3.1774 1 1])

%% Plot the residuals of L0
order = 0;

figure(2); clf
hold on; box on; grid on
[h,p] = residPlot(delta_store{order+1}, order, ind_col, colors);
ylim([0 .030])
xlim([-.6 .6])
set(gca,'fontsize',16)

%% Plot Order 2
order = 2;

figure(3); clf 
hold on; box on; grid on

[mdl, ~, ~] = special_plotter(b, delta_store, increment, ind_col, ind, colors, order);

% R^2 and P-value
text(8, .1, sprintf('R^2: %0.3f',mdl.Rsquared.Ordinary), 'fontsize', 12)
text(8, .07, sprintf('p-value: %0.3f',mdl.Coefficients.pValue(2)), 'fontsize', 12)

% Format
pbaspect([3.1774 1 1])
xlim([0 10])
ylim([-.3 .3])


%% Plot the residuals of L2 
order = 2;

figure(4); clf
hold on; box on
[h,p] = residPlot(delta_store{order+1}, order, ind_col, colors);
ylim([0 .03])
xlim([-.3 .3])

%% Plot Order 4
order = 4;

figure(5); clf 
hold on; box on; grid on

[mdl, ~, ~] = special_plotter(b, delta_store, increment, ind_col, ind, colors, order);

% R^2 and P-value
text(8, .1, sprintf('R^2: %0.3f',mdl.Rsquared.Ordinary), 'fontsize', 12)
text(8, .07, sprintf('p-value: %0.3f',mdl.Coefficients.pValue(2)), 'fontsize', 12)

% Format
pbaspect([3.1774 1 1])
xlim([0 10])
ylim([-.15 .15])


%% Some helper functions 
 function  [h,p] = residPlot(delta_store, order, ind_col, c)
    % Plot residuals
    normY= delta_store;
    
    for this_c = 2:length(ind_col)

        p_inds = [];

        % The inds of all the ooids from the sample
        if this_c == 2
            p_inds = ind_col(this_c-1):ind_col(this_c);
        else
            p_inds = ind_col(this_c-1)+1:ind_col(this_c);
        end  
        
        [h,p,~,~] = ks_test_normal(delta_store(p_inds))
        
        % Set bins
        lims = [min(delta_store(:)),max(delta_store(:))];
        range = lims(2)-lims(1);
        edges = linspace(lims(1)-0.1*range,lims(2)+0.1*range,30); % edges of bins
        centers = movmean(edges,2,'endpoints','discard'); % centers of bins

        gwin = gausswin(50); gwin = gwin/sum(gwin);
        hc = histcounts(delta_store(p_inds),edges); % Num in each bin
        xinterp = linspace(edges(1),edges(end),200); % make an x axis
        k = 5; % window for moving mean
        y = interp1(centers,movmean(hc,k,'endpoints','shrink'),xinterp); % interpolate
        y(isnan(y)) = 0; % remove any nans
        % plot
        %pY = fill(xinterp,conv(y/nansum(y),gwin,'same'),c(this_c,:),'facealpha',0,'edgecolor',c(this_c,:),'edgealpha',0.8);
        p = plot(xinterp, conv(y/nansum(y),gwin,'same'),'Color', c(this_c,:));
        p.LineWidth = 2;
    end
    plot([0 0], [0 3], 'k--')
    xlabel(sprintf('L_%i - L_%i', order, order))
    set(gca,'YTick', [], 'fontsize', 18)
    pbaspect([2 2 1])
end
 



function mdl = smallPlot(bs, ds)
    mdl = fitlm(bs(bs<2), ds(bs<2));
    p = plot(mdl); 
    delete(p(1)); delete(p(3)); delete(p(4)); delete(legend)
    p(2).Color = 'k';
    xlabel([]); ylabel([]); title([])
    pbaspect([2.1333, 1,1 ])
end
    
function [mdl, bs, ds, p] = special_plotter(b, delta_store, increment, ind_col, ind, colors, order)
    X = order+1;
    scat_scale = 100;
    bs = [];
    ds = [];

    for this_c = 2:length(ind_col)

        p_inds = [];

        % The inds of all the ooids from the sample
        if this_c == 2
            p_inds = ind_col(this_c-1):ind_col(this_c);
        else
            p_inds = ind_col(this_c-1)+1:ind_col(this_c);
        end  
        
        % The ooids that are too spherical
        p_inds_other = p_inds(~ind(p_inds));

        % the final indeces
        p_inds = p_inds(ind(p_inds));

        bs = [bs b{X}(p_inds)];
        ds  = [ds delta_store{X}(p_inds)];

        % plot the data that are not too spherical
        s1 = scatter(b{X}(p_inds), delta_store{X}(p_inds), scat_scale*increment{X}(p_inds), 'k', 'filled');
        s1.MarkerFaceColor = colors(this_c-1,:);
        s1.MarkerFaceAlpha = .8;

        % plot the very spherical data
        s = scatter(b{X}(p_inds_other),delta_store{X}(p_inds_other),scat_scale*increment{X}(p_inds_other));
        s.MarkerEdgeAlpha = .3;
        s.MarkerEdgeColor = colors(this_c-1,:);
        s.LineWidth = .0001;

    end

    % Plot the forecast line
    plot([0 100], [0 0], 'color', 'r', 'LineStyle', '--')

    mdl = fitlm(bs, ds);
    p = plot(mdl); 
    delete(p(1))
    delete(p(3))
    delete(p(4))
    delete(legend)
    p(2).Color = 'k';
    title([])
    ylabel(sprintf('L_%i - \hat{L}_%i', order, order), 'interpreter', 'tex')

    pbaspect([3.1774 1 1])
   

end



function mape = MAPE(Y,Y_pred)

    mape = (1/length(Y)) * sum(abs(Y - Y_pred)./abs(Y)) * 100;
    
end


function rsme = RMSE(Y,Y_pred)

%    mape = (1/length(Y)) * sum(abs(Y - Y_pred)./abs(Y)) * 100;
    rsme = sqrt(mean((Y - Y_pred).^2));

end



 function s = DDPlot(b, delta_store, y, increment, scat_scale, order)
    
    % Plot the data
    s = scatter(b, delta_store, scat_scale*increment, 'k', 'filled');
    s.MarkerFaceAlpha = .5;
    
    % Plot 0 line
    plot([0 max(b)+2], [0,0], 'r')

    % Set Labels
    xlabel('b [mm]')
    ylabel(sprintf('L_%i - \pred{L}_%i', order, order), 'interpreter', 'tex')

    keyboard
    if order == 0
        ylim([-1 1])
    end


 end

