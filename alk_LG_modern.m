% This script plots ooid occurences from the Recent and last glacial
% period. It plots their locations on a map of ocean alkalinity from the
% modern (from GLODAP Data) and from the last glacial (IPSL model)

clear; clc; close all;
%% Load data Geographic Datasets

% Ooid occurrences
[ooids_n] = readtable('data/ooid_occurrences.csv');
ooids_modern = ooids_n(ooids_n.age == 0,:);
ooids_nm = ooids_n(ooids_n.age == 0 & ooids_n.NormalMarine == 1,:); % modern, normal-marine ooids
ooids_nm = table2array(ooids_nm(:, [1 2]));

% Ooids from the last glacial period
ooids_gl = ooids_n(ooids_n.age >= 15 & ooids_n.age <= 120,:); % modern, normal-marine ooids
ooids_gl = table2array(ooids_gl(:, [1 2]));
%% Regions of Interest

% Global coral reef occurrences
% simplified dataset (created using the code in the lines below):
load('data/corals_simple.mat');

corals_bbox = nan(size(corals_simple,1),4); % holds lon1, lat1, lon2, lat2
for n = 1:size(corals_bbox,1)
    bbox = corals_simple(n).BoundingBox;
    corals_bbox(n,:) = [bbox(1,1) bbox(1,2), bbox(2,1), bbox(2,2)];
end
coral_locs = corals_bbox;

% GBB points:
lat_GBB = 23:0.2:27.5;
lon_GBB = -76:0.2:-73;
[TLAT,TLON] = meshgrid(lat_GBB,lon_GBB);
GBB_locs = [TLAT(:), TLON(:)];

% PG points:
lat_PG = 22:0.2:25;
lon_PG = 57:0.2:60;
[TLAT,TLON] = meshgrid(lat_PG,lon_PG);
PG_locs = [TLAT(:), TLON(:)];


% Eastern Mediterranean points:
lat_MD = 30:0.2:40;
lon_MD = 13:0.2:36;
[TLAT,TLON] = meshgrid(lat_MD,lon_MD);
MD_locs = [TLAT(:), TLON(:)];

% GLODAP Alkalinity dataset -- Modern
path = 'data/Alkalinity.nc';
[alk, GLODAP_lat, GLODAP_lon] = loadGLODAP(path, 'TAlk');
alk(alk < 0) = NaN;

% For whatever reason, the longitude are 20-380 in the GLODAP data, this resets them to -180
% 180. This regrids to match the GLODAP dataset
GLODAP_lon = GLODAP_lon-200;
alk_0m = repositionData(alk(:,:,1));
alk_100m = repositionData(alk(:,:,7));
alk_200m = repositionData(alk(:,:,10));

% It is gridded, so we inpaint for more complete coverage
alk_0m_inpaint = inpaintn(alk_0m);
alk_100m_inpaint = inpaintn(alk_100m);
alk_200m_inpaint = inpaintn(alk_200m);

%% Plot Modern Alkalinity Deviation from the mean
figure(2); clf
hold on
cmap = (brewermap(100,'RdYlBu'));

[LAT,LON] = setupMap(GLODAP_lat, GLODAP_lon);

% Calculate deviation from average value
var_plot = alk_200m - nanmean(alk_200m(:));

% The mediterranean is super alkaline, so it shouldnt set the limits of the
% colormap
crop_color = 220;
clm = [-1*max(abs(var_plot(:)))+crop_color max(abs(var_plot(:)))-crop_color];
formatMap(clm, '\deltaAlk (umol/kg)', cmap);
geoshow(LAT,LON, var_plot, 'displaytype', 'surface','alphadata',~isnan(var_plot));
plotm(ooids_nm(:,1),ooids_nm(:,2),6,'o','markerfacecolor','r','markeredgecolor','r');
gridm off
title('Dev. from modern mean ')

%% Load last glacial period alkalinity from IPSL model 

% Load
load('data/Alk_LGM')

% Getting the IPSL Model 
eval(['D=' 'I' ';'])
depth  = 200; 
[zc,pk]=(min(abs(D.z-depth(1)))); % Index of the correct depth

% modern - LGM
dAmodlgm=squeeze(diff(D.Alk(:,:,pk,:),[],4))';

% construct LGM by combining GLODAP and modeled LGM-Modern anomaly
try
    Algm=squeeze(GLO.Alk(:,:,pk))'-dAmodlgm; % LGM 
    Amod=squeeze(GLO.Alk(:,:,pk))';
    catch
end

% Transpose
Algm = Algm';

% Inpaint the lgm 
Algm_200m_inpaint = inpaintn(Algm);
LGM_dev_mod_bahamas = Algm - nanmean(alk_200m(:));


%% Plot Last glacial deviation from modern Bahamas
figure(4); clf
hold on

[LAT,LON] = setupMap(GLODAP_lat, GLODAP_lon);

var_plot = Algm_200m_inpaint- nanmean(alk_200m(:));
% The mediterranean is super alkaline, so it shouldnt set the limits of the
% colormap 
clm = [-1*max(abs(var_plot(:)))+crop_color max(abs(var_plot(:)))-crop_color];

formatMap(clm, '\deltaAlk (umol/kg)', cmap);
g = geoshow(LAT,LON-180, var_plot, 'displaytype', 'surface','alphadata',~isnan(Algm), 'FaceAlpha', 'interp');
plotm(ooids_gl(:,1),ooids_gl(:,2),6,'o','markerfacecolor','b','markeredgecolor','r');

title('LGM deviation from modern mean alk')

%% Sample parameters for locations

parameters_all = struct();
for n = 1:(size(coral_locs,1)+size(GBB_locs,1)+size(ooids_nm,1)+size(PG_locs,1)+size(MD_locs,1) + size(ooids_gl,1))
    % Draw values from entire ocean
    if n <= size(coral_locs,1)
        parameters_all(n).type = 0;
        parameters_all(n).lat = coral_locs(n,2);
        parameters_all(n).lon = coral_locs(n,1);
       
    % Draw values from coral locations
    elseif n <= size(coral_locs,1)+size(GBB_locs,1)
        parameters_all(n).type = 1;
        parameters_all(n).lat = GBB_locs(n-size(coral_locs,1),1);
        parameters_all(n).lon = GBB_locs(n-size(coral_locs,1),2);
        
    % Draw value from ooid locations
    elseif n <= size(coral_locs,1) + size(GBB_locs,1) + size(ooids_nm,1)
        parameters_all(n).type = 2;
        parameters_all(n).lat = ooids_nm(n-size(coral_locs,1)-size(GBB_locs,1),1);
        parameters_all(n).lon = ooids_nm(n-size(coral_locs,1)-size(GBB_locs,1),2);
        
    elseif n <= size(coral_locs,1) + size(GBB_locs,1) + size(ooids_nm,1) + size(PG_locs,1)
        parameters_all(n).type = 3;
        
        parameters_all(n).lat = PG_locs(n-size(coral_locs,1)-size(GBB_locs,1)-size(ooids_nm,1),1);
        parameters_all(n).lon = PG_locs(n-size(coral_locs,1)-size(GBB_locs,1)-size(ooids_nm,1),2);
        
    elseif n <= size(coral_locs,1) + size(GBB_locs,1) + size(ooids_nm,1) + size(PG_locs,1) +  size(MD_locs,1)
        
        parameters_all(n).type = 4;
        parameters_all(n).lat = MD_locs(n-size(coral_locs,1)-size(GBB_locs,1)-size(ooids_nm,1)-size(PG_locs,1),1);
        parameters_all(n).lon = MD_locs(n-size(coral_locs,1)-size(GBB_locs,1)-size(ooids_nm,1)-size(PG_locs,1),2);
    else
        parameters_all(n).type = 5;
        parameters_all(n).lat = ooids_gl(n-size(coral_locs,1)-size(GBB_locs,1)-size(ooids_nm,1)-size(PG_locs,1)-size(MD_locs,1),1);
        parameters_all(n).lon = ooids_gl(n-size(coral_locs,1)-size(GBB_locs,1)-size(ooids_nm,1)-size(PG_locs,1)-size(MD_locs,1),2);
   
    end
end
parameters_all = parameters_all';

for n = 1:length(parameters_all)
    % alkalinity:
    [~,indLat] = min(abs(GLODAP_lat - parameters_all(n).lat));
    [~,indLon] = min(abs(GLODAP_lon - parameters_all(n).lon));

    parameters_all(n).alk_0m = alk_0m(indLon,indLat);
    parameters_all(n).alk_100m = alk_100m(indLon,indLat);
    parameters_all(n).alk_200m = alk_200m(indLon,indLat); 
    parameters_all(n).LGM_alk_200m = Algm_200m_inpaint(indLon,indLat);

end

is_reef = [parameters_all(:).type]' == 0;
is_GBB = [parameters_all(:).type]' == 1;
is_ooid = [parameters_all(:).type]' == 2;
is_PG = [parameters_all(:).type]' == 3;
is_MD = [parameters_all(:).type]' == 4;
is_LGM_ooids = [parameters_all(:).type]' == 5;

%%
mask_am = 35;
lgm_polemask = Algm;
lgm_polemask(:,1:mask_am) = NaN;
lgm_polemask(:,180-mask_am:180) = NaN;


% Histograms
figure(6); clf
hold on; box off
c = (brewermap(8,'Set1'));
fields = fieldnames(parameters_all);
n = 6;
data_all = eval(fields{n});

% Mask out the polar oceans 
mod_polemask = data_all;
mod_polemask(:,1:mask_am) = NaN;
mod_polemask(:,180-mask_am:180) = NaN;

% set range
data_reefs = [parameters_all(is_reef).(fields{n})]';
lims = [min([data_all(:); data_reefs]),max([data_all(:); data_reefs])];
range = lims(2)-lims(1);
edges = linspace(lims(1)-0.1*range,lims(2)+0.1*range,100);
centers = movmean(edges,2,'endpoints','discard');

hc_all = histcounts(mod_polemask,edges);
hc_GBB = histcounts([parameters_all(is_GBB).(fields{n})]',edges);
hc_ooid = histcounts([parameters_all(is_ooid).(fields{n})]',edges);
hc_LGM_O = histcounts([parameters_all(is_LGM_ooids).LGM_alk_200m]',edges);
hc_pole_masked = 1;
hc_pole_masked_lgm = histcounts(lgm_polemask,edges);
hc_LGM_manual = histcounts([2352 2364 2366 2392 2394 2405 2407 2410 2489 2661],edges);

xinterp = linspace(edges(1),edges(end),200);

k = 5; % window for moving mean
y_all = interp1(centers,movmean(hc_all,k,'endpoints','shrink'),xinterp);
y_ooid = interp1(centers,movmean(hc_ooid,k,'endpoints','shrink'),xinterp);
y_LGM_O = interp1(centers,movmean(hc_LGM_O,k,'endpoints','shrink'),xinterp);
y_LGM_pm = interp1(centers,movmean(hc_pole_masked_lgm,k,'endpoints','shrink'),xinterp);
y_LGM_m = interp1(centers,movmean(hc_LGM_manual,k,'endpoints','shrink'),xinterp);

y_all(isnan(y_all)) = 0;
y_ooid(isnan(y_ooid)) = 0;
y_LGM_O(isnan(y_LGM_O)) = 0;
y_LGM_pm(isnan(y_LGM_pm)) = 0;
y_LGM_m(isnan(y_LGM_m)) = 0;

% smoothing 
gwin = gausswin(10); gwin = gwin/sum(gwin);

% plotting colors
colA = hex2rgb('#fa9442');
colLGM = hex2rgb('#008aa1');

% Modern Ocean
pA = fill(xinterp,conv(y_all/nansum(y_all),gwin,'same'),colA,'facealpha',0,'edgecolor',colA, 'LineWidth', 4);

% Plot polemasked 
LGM_pm = fill(xinterp,conv(y_LGM_pm/nansum(y_LGM_pm),gwin,'same'),colLGM,'facealpha',0,'edgecolor',colLGM, 'LineWidth',4);

% modern ooids
mooids = max([nansum(y_ooid) nansum(y_LGM_O)]);
pO = fill(xinterp,conv(y_ooid/nansum(y_ooid),gwin,'same'),colA,'facealpha',0.2,'edgecolor',colA, 'LineWidth',2);

% LGM ooids
gwin = gausswin(10); gwin = gwin/sum(gwin);
LGM_po = fill(xinterp,conv(y_LGM_m/nansum(y_LGM_m),gwin,'same'),colLGM,'facealpha',0.2,'edgecolor',colLGM, 'LineWidth',2);

% Format
xlim([2200 2550]);
set(gca,'ytick',[]);
xlabel('Alk (umol/kg)')
set(gca, 'Fontsize', 14)
pbaspect([.95 .75 1]) 


%% Helper functions

function dataset_final = repositionData(dataset)

    % Data are offset  by 160 degrees in longitude. This repositions it
    for x = 1:360
        n = (x-160);
        if n < 1
            n = n+360;
        end
        dataset_final(n,:) = dataset(x,:);

    end
end


function [LAT, LON] = setupMap(GLODAP_lat, GLODAP_lon)

    load coastlines % quick load in coastlines
    worldmap('World')
    geoshow(coastlat,coastlon,'DisplayType','polygon','FaceColor',[1 1 1]/1.2,'edgecolor','k');
    [LAT,LON] = meshgrid(GLODAP_lat, GLODAP_lon); % Plotting grid for GLodap

end

function formatMap(clims, ylab, cmap)

                     
    cb = colorbar('southoutside');      % Put colorbar on the bottom
    colormap(cmap);
    ylabel(cb,ylab);                    % Label colorbar
    caxis([clims]);                     % Set color limits
    gridm off; mlabel off;              % removes the longitude labels
    set(gca, 'fontsize', 16)            % Set fontsize
end