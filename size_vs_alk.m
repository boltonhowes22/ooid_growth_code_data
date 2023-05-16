clear; clc; close all;
save_path = '../';
c = [0 0 0;
    0 0 0;
    0 0 0];

%% Load data Geographic Datasets

% 0. Matlab global coastlines:
load coastlines

% 1. Ooid occurrences
ooids_n = readtable('data/ooid_occurrences.csv');
ooids_modern = ooids_n(ooids_n.age == 0,:);
ooids_nm = ooids_n(ooids_n.age == 0 & ooids_n.NormalMarine == 1,:); % modern, normal-marine ooids
ooids_gl = ooids_n(ooids_n.age >= 15 & ooids_n.age <= 120,:); % modern, normal-marine ooids

% 2. Global coral reef occurrences
% simplified dataset (created using the code in the lines below):
load('data/corals_simple.mat');

corals_bbox = nan(size(corals_simple,1),4); % holds lon1, lat1, lon2, lat2
for n = 1:size(corals_bbox,1)
    bbox = corals_simple(n).BoundingBox;
    corals_bbox(n,:) = [bbox(1,1) bbox(1,2), bbox(2,1), bbox(2,2)];
end

%% Regions of Interest
coral_locs = corals_bbox;

% GBB points:
lat_GBB = 23:0.2:27.5;
lon_GBB = -76:0.2:-73;
[TLAT,TLON] = meshgrid(lat_GBB,lon_GBB);
GBB_locs = [TLAT(:), TLON(:)];
%GBB_locs = GBB_locs';

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

% Reunion
lat_R = -22:0.2:-18;
lon_R = 63:0.2:65;
[TLAT,TLON] = meshgrid(lat_R,lon_R);
R_locs = [TLAT(:), TLON(:)];

% Bengal Shelf
lat_B = 15:0.2:19;
lon_B = 84:0.2:88;
[TLAT,TLON] = meshgrid(lat_B,lon_B);
B_locs = [TLAT(:), TLON(:)];

% Georgia
lat_G = 30:0.2:32;
lon_G = -81:0.2:-79;
[TLAT,TLON] = meshgrid(lat_G,lon_G);
G_locs = [TLAT(:), TLON(:)];

% NW India
lat_NW = 19.5:0.2:22.5;
lon_NW = 66:0.2:71;
[TLAT,TLON] = meshgrid(lat_NW,lon_NW);
NW_locs = [TLAT(:), TLON(:)];

% Shark Bay
lat_SB = -26:0.2:23;
lon_SB = 109:0.2:113;
[TLAT,TLON] = meshgrid(lat_SB,lon_SB);
SB_locs = [TLAT(:), TLON(:)];

%% 3. GLODAP Alkalinity dataset
path = 'data/Alkalinity.nc';
[alk, GLODAP_lat, GLODAP_lon] = loadGLODAP(path, 'TAlk');
alk(alk < 0) = NaN;

% For whatever reason, the longitude are 20-380, this resets them to -180
% 180
GLODAP_lon = GLODAP_lon-200;

alk_0m = repositionData(alk(:,:,1));
alk_100m = repositionData(alk(:,:,7));
alk_200m = repositionData(alk(:,:,10));

alk_0m_inpaint = inpaintn(alk_0m);
alk_100m_inpaint = inpaintn(alk_100m);
alk_200m_inpaint = inpaintn(alk_200m);


%% Plot LGM Dev from modern Bahamas
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
%

Algm = Algm';
Algm_200m_inpaint = inpaintn(Algm);


%% Sample parameters for locations
sz_cl = size(coral_locs,1);
sz_gbb = sz_cl + size(GBB_locs,1);  
sz_nm = sz_gbb +size(ooids_nm,1);
sz_pg = sz_nm + size(PG_locs,1);
sz_MD = sz_pg +size(MD_locs,1);
sz_gl = sz_MD +  size(ooids_gl,1);
sz_R = sz_gl + size(R_locs, 1);
sz_B = sz_R + size(B_locs,1);
sz_G =  sz_B + size(G_locs,1);
sz_NW = sz_G + size(NW_locs,1);
sz_SB = sz_NW + size(SB_locs,1);

parameters_all = struct();
for n = 1:sz_SB%(size(coral_locs,1)+size(GBB_locs,1)+size(ooids_nm,1)+size(PG_locs,1)+size(MD_locs,1) + size(ooids_gl,1))
    % Draw values from entire ocean
    if n <= sz_cl
        parameters_all(n).type = 0;
        parameters_all(n).lat = coral_locs(n,2);
        parameters_all(n).lon = coral_locs(n,1);
        
    % Draw values from coral locations
    elseif n <= sz_gbb
        parameters_all(n).type = 1;
        parameters_all(n).lat = GBB_locs(n-sz_cl,1);
        parameters_all(n).lon = GBB_locs(n-sz_cl,2);
        
    % Draw value from ooid locations
    elseif n <= sz_nm
        
        parameters_all(n).type = 2;
        parameters_all(n).lat = ooids_nm(n-sz_gbb,1).latitude;
        parameters_all(n).lon = ooids_nm(n-sz_gbb,2).longitude;
        
    elseif n <= sz_pg  
        
        parameters_all(n).type = 3;
        parameters_all(n).lat = PG_locs(n-sz_nm,1);
        parameters_all(n).lon = PG_locs(n-sz_nm,2);
        
    elseif n <= sz_MD 
        
        parameters_all(n).type = 4;
        parameters_all(n).lat = MD_locs(n-sz_pg,1);
        parameters_all(n).lon = MD_locs(n-sz_pg,2);
    
    elseif n <= sz_gl 
    
        parameters_all(n).type = 5;
        parameters_all(n).lat = ooids_gl(n-sz_MD,1).latitude;
        parameters_all(n).lon = ooids_gl(n-sz_MD,2).longitude;
        
    elseif n <= sz_R % now Rodrigues
        
        parameters_all(n).type = 6;
        parameters_all(n).lat = R_locs(n-sz_gl,1);
        parameters_all(n).lon = R_locs(n-sz_gl,2);
        
    elseif n <= sz_B % now Bengal Shelf
        
        parameters_all(n).type = 7;
        parameters_all(n).lat = B_locs(n-sz_R,1);
        parameters_all(n).lon = B_locs(n-sz_R,2);
        
    elseif n <= sz_G % now Georgia
        
        parameters_all(n).type = 8;
        parameters_all(n).lat = G_locs(n-sz_B,1);
        parameters_all(n).lon = G_locs(n-sz_B,2);
        
    elseif n <= sz_NW % now NW india
        
        parameters_all(n).type = 9;
        parameters_all(n).lat = NW_locs(n-sz_G,1);
        parameters_all(n).lon = NW_locs(n-sz_G,2);

    elseif n <= sz_SB % now NW india
        
        parameters_all(n).type = 10;
        parameters_all(n).lat = SB_locs(n-sz_NW,1);
        parameters_all(n).lon = SB_locs(n-sz_NW,2);
              
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
    parameters_all(n).lgm_alk_200m = Algm_200m_inpaint(indLon,indLat);
    

end

is_reef = [parameters_all(:).type]' == 0;
is_GBB  = [parameters_all(:).type]' == 1;
is_ooid = [parameters_all(:).type]' == 2;
is_PG   = [parameters_all(:).type]' == 3;
is_MD   = [parameters_all(:).type]' == 4;
is_LGM_ooids = [parameters_all(:).type]' == 5;
is_R    = [parameters_all(:).type]' == 6;
is_B    = [parameters_all(:).type]' == 7;
is_G    = [parameters_all(:).type]' == 8;
is_NW   = [parameters_all(:).type]' == 9;
is_SB   = [parameters_all(:).type]' == 10;

%% Plot ooid size distributions as a function of alkalinity 

% andros = 2375;
em_andros = [350 500] ; % microns get exact numbers
pers_gulf = [500 390];
tc = [430 475];
shark_bay = [253, 365];
alexandria = [413, 535];
bengal_shelf = [300 400 350 450 500];
georgia = [330 380 450];
nw_india = [700 900];

% modern Data points
data_GBB = [parameters_all(is_GBB).('alk_200m')]';
data_MD = [parameters_all(is_MD).('alk_200m')]';
data_PG = [parameters_all(is_PG).('alk_200m')]';
data_SB = [parameters_all(is_SB).('alk_200m')]';

% LGM Data points
data_B = [parameters_all(is_B).('lgm_alk_200m')]';
data_G = [parameters_all(is_G).('lgm_alk_200m')]';


f = figure(14);clf
hold on; box on; grid on;
c = brewermap([6],'Accent');
ylim([250 550])

% Plot the size as a function of alkalinity
ms = 50;
sizeAlkPlot(data_GBB, em_andros, 'k')
sizeAlkPlot(data_PG, pers_gulf, 'k')
sizeAlkPlot(data_GBB, tc, 'k')
sizeAlkPlot(data_MD, alexandria, 'k')
sizeAlkPlot(data_SB, shark_bay, 'k')
sizeAlkPlot(data_B, bengal_shelf, [1 1 1]./1.7)
sizeAlkPlot(data_G, georgia, [1 1 1]./1.7)


% Plot the global ocean and the coral reef alkalinity 
data_all = eval('alk_200m');
plot([nanmean(data_all(:)) nanmean(data_all(:))], [0 10000], 'k')
reef_all = [parameters_all(is_reef).('alk_200m')];
plot([nanmean(reef_all(:)) nanmean(reef_all(:))], [0 10000], 'k')
text(nanmean(reef_all(:))+8, 475, 'Coral Reef Mean', 'rot', 90)
text(nanmean(data_all(:))-8, 475, 'Global Mean', 'rot', 90)


% Linear regression
all_size  = [mean(shark_bay), mean(bengal_shelf), mean(georgia),...
            mean(em_andros), mean(alexandria), mean(pers_gulf), mean(tc)];
all_alk =   [nanmean(data_SB), nanmean(data_B), nanmean(data_G),...
            nanmean(data_GBB), nanmean(data_MD), nanmean(data_PG), nanmean(data_GBB)];
 
[X, s_ind] = sort(all_alk);
y = all_size(s_ind);
p = polyfit(log(X), log(y), 1);
Xf = min(X):max(X);
y = polyval(p, log(Xf));

% Plot the line of best fit 
plot(Xf, exp(y), 'k')

% Format
ylabel('Diameter (um)', 'fontsize', 14)
xlabel('Alkalinity (umol/kg)', 'fontsize', 14)
pbaspect([1 1 1])


%% Let's see what percentiles the ooids are in for the entire ocean 

fields = fieldnames(parameters_all);
% Print the percentiles for the localities
for i = 4:length(fields)-1
    
    data_reefs = [parameters_all(is_reef).(fields{i})]';
    data_ooids = [parameters_all(is_ooid).(fields{i})]';
    reef_area = [corals_simple(:).Shape_Area]';
    
    data_all = eval(fields{i});
    data_GBB = [parameters_all(is_GBB).(fields{i})]';

    % where does the mean of the Bahamas show up in the global dataset:
    val = nanmean(data_ooids);
    prct = 100*sum((data_reefs < val).*reef_area/sum(reef_area)); % weighted by reef area
    
    fprintf('%s (%s = %2.2f): %2.2fth percentile\n',fields{i},'\mu',val,prct);
end


%% Helper Functions
function sizeAlkPlot(alk_data, size_data, color)
    ms = 125;
    scatter(nanmean(alk_data), mean(size_data), ms, color(1,:), 'filled')
    plot([nanmean(alk_data) nanmean(alk_data)], [prctile(size_data,25) prctile(size_data,75)],...
        'Color', color(1,:), 'LineWidth', 2.5)
  

end


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



