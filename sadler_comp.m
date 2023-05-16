colors = [hex2rgb('#084887');
          hex2rgb('#8B575C');
          hex2rgb('#4B644A');
          hex2rgb('#B2B1CF');
          hex2rgb('#C1EEFF');];
msize = 50; % markersize
lw = 2;

density = 2.93; % Density of Aragonite


highborne_sz = .0250; %cm diameter
highborne_ls = 800; %yrs
h1 = growthRate(highborne_sz, highborne_ls);

carbla_sz1 = .0250; %cm diameter
carbla_ls1 = 1470; %yrs
c1 = growthRate(carbla_sz1, carbla_ls1);


carbla_sz2 = .0590; %cm diameter
carbla_ls2 = 1120; %yrs
c2 = growthRate(carbla_sz1, carbla_ls1);

%% Ferguson 1978 Rates: Synthetic ooid growth
% experiments lasted between 1 an 30 days. we chose a middle value of 15.

% They were able to grow synthetic ooids of at least 400 microns in
% diameter
% The minumum rate is for if it took 30 days. The maximum rate is for
% taking 1 day. 
diameter = .03; % cm
ferguson_timescale = [15]./365; % years
ferguson_rate = growthRate(diameter, ferguson_timescale);

%% The Broecker and Takahashi Data
% precipitation rate 
bro_taka_rate = [140 44 32  68 62 22 58 53 82 70 75 85 38 46 62 42];

% Age of the water mass
bro_taka_age_years =  [20 150 245 12 6 22 48 71 56 73 66 56 142 133 96 180]./365; % age in years
BT = BroeckerConverter(bro_taka_rate, density);

%% Rate from Smith 
% Canton Atoll 
% .5 kg/m^2/yr over 50 days on average
canton_rate = SmithConverter(0.5, density);
canton_timespan = 50; % days, reported mean. Maximum is 95 days, based on salinity 
canton_timespan = canton_timespan/365; % in years


%% Davies and Kinsey 1977
% rates in kg/m2/yr
algal_flat = 4; 
reef_flat = 4.5;
sand_flast = .3;
retic_path_reef_lagoon = 1.5;
deep_lagoon = 0.5;
integrated = 1.5;
OTR_envs = [algal_flat,reef_flat, sand_flast, retic_path_reef_lagoon, deep_lagoon, integrated];

for x = 1:length(OTR_envs)
    OTR_rates(x) = SmithConverter(OTR_envs(x), density); % Same conversion as Smith's data
end

% values is the average annual rate monitored each year for 6 years
OTR_timespan = 1; %year

%% Bialik 2022 Mediterranean
% The 459 day and 50 measurement are the runtime for the entire experiment 
bialik_timespan = [33 33 43 29 459]./365;
bialik_measurements = [60 120 250 220 50];

for x = 1:length(bialik_measurements)  
    bialik_rates(x) = BialikConverter(bialik_measurements(x), density);
end

%% Diurnal engine peak rate
% using parameters from Geyman, Diurnal Engine 

k = log(1.11); %micromol/m2/hr
OmA = 6.5; % Afternoon saturation state
n = 2.2; 
F = k*(OmA -1)^n; % micromol/m2/hr
molarmass = 100.087;

% Convert to moles
F_moles = F * (1*10^-6);


duration = 3; % Let's say this  continues for 3 hours in the afternoon
converted = (((24/duration)*365*F) * (1*10^-6) * molarmass) / density;
afternoon = duration/24/365; % 3 hours

% Plotting
figure(1); clf
hold on; box on; grid off;
xlabel('Measurement Interval (yrs)')
ylabel('Precipitation Rate (cm/yr)')

cm = brewermap([7],'Accent');

msize1 = 70;
alpha = .75;
% Plot the Ooids 
% radiocarbon rates
sc = scatter(highborne_ls, h1, msize1, 'o');
sc.LineWidth = 2;
sc = scatter([carbla_ls1 , carbla_ls2], [c1, c2], msize1, 'o');
sc.LineWidth = 2;

% Plot the Broecker and Takahashi Rates
sc = scatter(bro_taka_age_years, BT, msize1, 'filled');
sc.MarkerFaceColor = cm(1,:);
sc.MarkerFaceAlpha = alpha;

% Plot the Smith Rates from Canton Atoll
sc = scatter(canton_timespan, canton_rate, msize1,'filled');
sc.MarkerFaceColor = cm(2,:);
sc.MarkerFaceAlpha = alpha;

% Plot Davies and Kinsey Rates from One Tree Reef
sc = scatter(repmat(OTR_timespan,[1 length(OTR_envs)]), OTR_rates, msize1, 'filled');
sc.MarkerFaceColor = cm(3,:);
sc.MarkerFaceAlpha = alpha;

% Plot Ferguson rates
sc = scatter(ferguson_timescale, ferguson_rate, msize1, 'filled');
sc.MarkerFaceColor = cm(7,:);
sc.MarkerFaceAlpha = alpha;

% Plot Bialik rates
sc = scatter(bialik_timespan, bialik_rates,msize1, 'filled');
sc.MarkerFaceColor = cm(5,:);
sc.MarkerFaceAlpha = alpha;

% Plot diurnal engine estimate
sc = scatter(afternoon, converted, msize1, 'filled');
sc.MarkerFaceColor = cm(6,:);
sc.MarkerFaceAlpha = alpha;

% Load in speleothem data
path_to_sisalv2 = 'data/dating.csv';
[spel_age_diff, spel_growth_rate] = speleothemDatasetSadler(path_to_sisalv2);

% Plot subset of speleothem data 
plot_i = randi(length(spel_growth_rate), [1, 50000]);
sc = scatter(spel_age_diff(plot_i), spel_growth_rate(plot_i), 'filled', 'k');
sc.SizeData = 10;
sc.MarkerFaceAlpha = 1;
sc.MarkerEdgeColor = 'None';

% Plot the Fit to the data
% Compile data for line of best fit
[time_scale, ind] = sort([bialik_rates,highborne_ls, carbla_ls1, carbla_ls2, bro_taka_age_years, canton_timespan, repmat(OTR_timespan,[1 length(OTR_envs)]), afternoon]);
Y = [bialik_timespan,h1, c1, c2, BT, canton_rate, OTR_rates, converted];
Y_sorted = Y(ind);

% Fit in log-log space
% assist from https://www.mathworks.com/matlabcentral/answers/1439604-best-fit-line-in-log-log-scale
p = polyfit(log(time_scale), log(Y_sorted), 1);
y = polyval(p, log(time_scale));

% Plot the line of best fit 
plot(time_scale, exp(y), 'k')

[time_scale, ind] = sort(spel_age_diff);
Y_sorted = spel_growth_rate(ind);

%scatter(log(all_mat(:,1)), log(all_mat(:,2)))
% Fit in log-log space
% assist from https://www.mathworks.com/matlabcentral/answers/1439604-best-fit-line-in-log-log-scale
p = polyfit(log(spel_age_diff), log(spel_growth_rate), 1);
y = polyval(p, log(spel_age_diff(:,1)));

% Plot the line of best fit 
plot(spel_age_diff(:,1), exp(y), 'k')



gridx = [3/365 1  100 10000];
ts_text = {'3 day'; '1 yr'; '100 yrs'; '10 kyr'};
for x = 1:length(gridx)
   p = plot([gridx(x) gridx(x)], [10^-20 10^8]); 
   p.LineWidth = .1;
   p.Color = [.1 .1 .1 .8];
   t = text(gridx(x)-gridx(x)*.3, .15e-5,ts_text{x}, 'Rotation', 90, 'fontsize', 12);
end

xlim([10^-4 15^5])
ylim([12^-6 8^1])
set(gca,'YScale', 'log', 'XScale', 'log', 'Fontsize', 14)
pbaspect([3 3 1])

%savePlot('/Users/boltonhowes/Dropbox (Princeton)/ooids/ooid_growth_repo/LaTeXTemplates_pnas_v1.44/manuscript_figures/discussion_fig/', 'rate_sadler_sisalv.eps')
%%
function cm3cm2yr = BroeckerConverter(mgcm2yr, density)
% Convert mg/cm^2/yr to cm^3/cm^2/yr
% Inputs:
% mgcm2yr: mass deposition rate in mg/cm^2/yr
% density: density of the material in g/cm^3
% Output:
% volume: volume deposition rate in cm^3/cm^2/yr

% Calculate volume deposition rate in cm^3/cm^2/yr
cm3cm2yr = mgcm2yr * 0.0001 / density;

end


function cm3cm2yr = SmithConverter(kgm2yr, density)
% Convert kg/m^2/yr to cm^3/cm^2/yr
% Inputs:
% kgm2yr: mass deposition rate in kg/m^2/yr
% density: density of the material in g/cm^3
% Output:
% volume: volume deposition rate in cm^3/cm^2/yr

% Calculate volume deposition rate in cm^3/cm^2/yr
cm3cm2yr = (kgm2yr * 0.01) / density ;

end


function growth_rate = growthRate(diameter, lifespan)
% Calculate growth rate in cm/yr
%
% In
% diameter: in cm 
% lifespan: in years
% 
% out
% growth_rate: in micron^3/hr
% Begin %%%%

    radius = (diameter/2);    
    growth_rate = radius/lifespan;
end


function cm3cm2yr = BialikConverter(mgm2day, density)

    % convert from per day to year
    mgm2year = mgm2day * 365;

    % convert from per m2 to per cm2
    mgcm2year = mgm2year/10000;
   
    % convert from mg to g
    gcm2year = mgcm2year/1000;

    % convert from g to cm3
    cm3cm2yr = gcm2year/density;    
end


function [s, all_mat] = plotSpeleothem(path_name, name)

    % load in and plot the speleothem data
    dime3 = readtable([path_name name], 'Sheet', 'dime3');
    dime2 = readtable([path_name name], 'Sheet', 'dime2');
    dime4 = readtable([path_name name], 'Sheet', 'dime4');
    h12 = readtable([path_name name], 'Sheet', '7H12');
    xt4 = readtable([path_name name], 'Sheet', 'XT4');
    xt5 = readtable([path_name name], 'Sheet', 'XT5');
    BF2_1 = readtable([path_name name], 'Sheet', 'BF2_1');
    BF2_2 = readtable([path_name name], 'Sheet', 'BF2_2');

    all = {dime3; dime2; dime4; h12; xt4; xt5; BF2_1; BF2_2};
    
    for x = 1:length(all)
       s{x} = simpleDist(all{x});

    end
   
    all_mat = [];
   
    for d = 1:length(all)
        
        all_mat = [all_mat; s{d}];
    end
    
end


function s = simpleDist(df)
    
    % Just calculate the growth rate for the speleothem data
    distance = df.DFT_mm * .1; % convert to cm 
    age = df.age_kyr*1000; % conver to years
    %age = age*365*24; % convert to hours
    age_diff = diff(age);
    gr = diff(distance)./diff(age*365*24);
    
    ind = age_diff > 0;
    s = [age_diff(ind) gr(ind)];
end




