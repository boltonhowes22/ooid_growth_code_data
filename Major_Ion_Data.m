% Emily Geyman
% Created: 01.06.23

% Make a Phanerozoic curve of Mg, Ca based on the data:

clear; clc; close all;

%% Load data
% MgCa:
[dn,ds,~] = xlsread('Datasets/Major_Ion_Concentrations.xlsx','Rees2010_Phanerozoic');
isFi = strcmp(ds(2:end,5),'Fluid_Inclusion');
isEc = strcmp(ds(2:end,5),'Echinoid');
MgCa_Fi = dn(isFi,1:4); % age, low, med, high
MgCa_Ec = dn(isEc,1:4); % age, low, med, high
[MgCa,ds,~] = xlsread('Datasets/Major_Ion_Concentrations.xlsx','MgCa');
isD = strcmp(ds(2:end,3),'Dickson2002');
isC = strcmp(ds(2:end,3),'Coggon2010');
isR = strcmp(ds(2:end,3),'Rausch2013');
isG = strcmp(ds(2:end,3),'Gothman2015');


figure(1)
clf
subplot(3,1,1);
hold on; box on; grid on;
% Fluid inclusion:
for i = 1:size(MgCa_Fi,1)
    plot([MgCa_Fi(i,1),MgCa_Fi(i,1)],[MgCa_Fi(i,2) MgCa_Fi(i,4)],'-','linewidth',2,'color',[1 1 1]/2);
end
pFi = plot(MgCa_Fi(:,1),MgCa_Fi(:,3),'ks','markerfacecolor',[1 1 1]/1.5,'markersize',8);
% Echinoid
for i = 1:size(MgCa_Ec,1)
    plot([MgCa_Ec(i,1),MgCa_Ec(i,1)],[MgCa_Ec(i,2) MgCa_Ec(i,4)],'-','linewidth',2,'color',[1 1 1]/1.2);
end
pFi = plot(MgCa_Ec(:,1),MgCa_Ec(:,3),'ks','markerfacecolor',[1 1 1]/1.2,'markersize',8);
pD = plot(MgCa(isD,1),MgCa(isD,2),'k^','markerfacecolor',[1 1 1]/1.3,'markersize',8);
pC = plot(MgCa(isC,1),MgCa(isC,2),'kd','markerfacecolor',[1 1 1]/1.3,'markersize',8);
pR = plot(MgCa(isR,1),MgCa(isR,2),'kp','markerfacecolor',[1 1 1]/1.3,'markersize',8);
%pG = plot(MgCa(isG,1),MgCa(isG,2),'kv','markerfacecolor',[1 1 1]/1.3,'markersize',8);
pM = plot(0,5.09,'kp','markersize',14,'markerfacecolor',[1 1 1]/4);
set(gca,'xdir','reverse');
ylabel('Mg/Ca');

[Holt_Ca,~,~] = xlsread('Datasets/Major_Ion_Concentrations.xlsx','Holt2014_Phanerozoic_Ca');
subplot(3,1,2);
hold on; box on; grid on;
for i = 1:size(Holt_Ca,1)
    plot([Holt_Ca(i,1),Holt_Ca(i,1)],[Holt_Ca(i,2),Holt_Ca(i,3)],'-','linewidth',4,'color',[1 1 1]/1.8);
end
pM = plot(0,10.2,'kp','markersize',14,'markerfacecolor',[1 1 1]/4);
set(gca,'xdir','reverse');
ylabel('[Ca] [mmol/kg]');

[Antonelli_Mg,~,~] = xlsread('Datasets/Major_Ion_Concentrations.xlsx','Mg_Phanerozoic_Antonelli');
subplot(3,1,3);
hold on; box on; grid on;
for i = 1:size(Antonelli_Mg,1)
    plot([Antonelli_Mg(i,1),Antonelli_Mg(i,1)],[Antonelli_Mg(i,2),Antonelli_Mg(i,4)],'-','linewidth',4,'color',[1 1 1]/1.8);
end
pM = plot(0,52,'kp','markersize',14,'markerfacecolor',[1 1 1]/4);
set(gca,'xdir','reverse');
ylabel('[Mg] [mmol/kg]');

%% Do an MCMC approach to a joint inversion of the data

% Prepare organized datasets:
MgCa_all = [MgCa_Fi; MgCa_Ec; [MgCa(~isG,1),NaN*MgCa(~isG,2),MgCa(~isG,2),NaN*MgCa(~isG,2)]];
% assign average erorr bars to the data points w/o errors:
% error bar size as a function of Mg/Ca:
er = (MgCa_all(:,4)-MgCa_all(:,2)) / 2;
v = ~isnan(er);
X = polyfit(MgCa_all(v,3),er(v),1);
er_fillmissing = polyval(X,MgCa_all(~v,3));
MgCa_all(~v,[2,4]) = [MgCa_all(~v,3)-er_fillmissing,MgCa_all(~v,3)+er_fillmissing];

% remove outliers:
outliers = [8.52,5.46; 23.65,3.63; 100.13,0.65; 130.74,1.85; 220.49,4.55; ...
    384.44,2.56; 515.54,1.18; 250.11,3.56; 253.36,3.74; 302.72,2.89];
thresh_outlier = 0.05;
d = min((MgCa_all(:,1)-outliers(:,1)').^2 + (MgCa_all(:,3)-outliers(:,2)').^2,[],2);
is_outlier = d < thresh_outlier;

MgCa_all = MgCa_all(~is_outlier,:);

%%

prc_intervals = 10:5:90;
x = norminv(prc_intervals/100);
%cmap = flipud(brewermap(2*floor(length(x)/2),'BrBG'));
%cmap = flipud(cmap(1:floor(length(x)/2),:));
cmap = brewermap(floor(length(x)/2)+1,'Blues');

xlims = [0 570];
xpred = linspace(600,0,200)';

figure(2)
clf
subplot(3,1,1);
hold on; box on; grid on;
set(gca,'xdir','reverse');
ylabel('Mg/Ca');
% fit the Mg/Ca data:
gprMdl = fitrgp(MgCa_all(:,1),MgCa_all(:,3),'KernelFunction','exponential');
[ypred,ysd,yint] = predict(gprMdl,xpred);
%fill([xpred; flipud(xpred)],[ypred-ysd; flipud(ypred+ysd)],[1 1 1]/3,'facealpha',0.2);
for i = 1:floor(length(x)/2)
    fill([xpred; flipud(xpred)],[ypred-abs(x(i))*ysd; flipud(ypred+abs(x(i))*ysd)],cmap(i,:),'edgecolor','none');
end
%plot(xpred,ypred,'k-','linewidth',2);
%for i = 1:size(MgCa_all,1)
    %plot([MgCa_all(i,1),MgCa_all(i,1)],[MgCa_all(i,2),MgCa_all(i,4)],'k-');
%end
%plot(MgCa_all(:,1),MgCa_all(:,3),'ks','markerfacecolor',[1 1 1]/3);
pM = plot(0,5.09,'kp','markersize',14,'markerfacecolor',[1 1 1]/3);
xlim(xlims);
ylim([0 6]);

% what happens if you move Ca and Mg 1:1 from starting values following
% Mg/Ca fit?
num_iter = 100;
Ca_fit1 = nan(length(ypred),num_iter);
Mg_fit1 = nan(length(ypred),num_iter);
Ca_fit1(end,:) = 10.2;
Mg_fit1(end,:) = 52;
rng('shuffle');
for ii = 1:num_iter
    offset_iter = randn(1);
    for i = 1:length(ypred)-1
        ratio_pick = ypred(end-i) + ysd(end-i)*offset_iter;
        delta = fminsearch(@(delta)abs((Mg_fit1(end-(i-1),ii)-delta)/(Ca_fit1(end-(i-1),ii)+delta)-ratio_pick),0);
        Ca_fit1(end-i,ii) = Ca_fit1(end-(i-1),ii)+delta;
        Mg_fit1(end-i,ii) = Mg_fit1(end-(i-1),ii)-delta;
    end
end
for i = 1:size(MgCa_Fi,1)
    plot([MgCa_Fi(i,1),MgCa_Fi(i,1)],[MgCa_Fi(i,2) MgCa_Fi(i,4)],'-','linewidth',1,'color',[1 1 1]/2);
end
pFi = plot(MgCa_Fi(:,1),MgCa_Fi(:,3),'ks','markerfacecolor',[1 1 1]/1.5,'markersize',5);
% Echinoid
for i = 1:size(MgCa_Ec,1)
    plot([MgCa_Ec(i,1),MgCa_Ec(i,1)],[MgCa_Ec(i,2) MgCa_Ec(i,4)],'-','linewidth',1,'color',[1 1 1]/1.2);
end
pEc = plot(MgCa_Ec(:,1),MgCa_Ec(:,3),'ko','markerfacecolor',[1 1 1]/1.2,'markersize',5);
pD = plot(MgCa(isD,1),MgCa(isD,2),'k^','markerfacecolor',[1 1 1]/1.3,'markersize',5);
pC = plot(MgCa(isC,1),MgCa(isC,2),'kd','markerfacecolor',[1 1 1]/1.3,'markersize',5);
pR = plot(MgCa(isR,1),MgCa(isR,2),'kp','markerfacecolor',[1 1 1]/1.3,'markersize',5);
% legend([pFi pEc pD pC pR],'Ries (2010): Fluid inclusions','Ries (2010): Echnoids',...
%     'Dickson (2002)','Coggon et al. (2010)','Rausch et al. (2013)','location','eastoutside');
% cb = colorbar;
% colormap([cmap(1:end-1,:); flipud(cmap(1:end-2,:))]);

subplot(3,1,2);
hold on; box on; grid on;
for i = 1:floor(length(x)/2)
    fill([xpred; flipud(xpred)],[prctile(Ca_fit1,prc_intervals(i),2); flipud(prctile(Ca_fit1,prc_intervals(end+1-i),2))],cmap(i,:),'edgecolor','none');
end
for i = 1:size(Holt_Ca,1)
    plot([Holt_Ca(i,1),Holt_Ca(i,1)],[Holt_Ca(i,2),Holt_Ca(i,3)],'-','linewidth',4,'color',[1 1 1]/3);
end
pM = plot(0,10.2,'kp','markersize',14,'markerfacecolor',[1 1 1]/3);
% for i = 1:num_iter
%     plot(xpred,Ca_fit1(:,i),'k-');
% end
set(gca,'xdir','reverse');
ylabel('[Ca] [mmol/kg]');
xlim(xlims);
ylim([0 60]);

subplot(3,1,3);
hold on; box on; grid on;
for i = 1:floor(length(x)/2)
    fill([xpred; flipud(xpred)],[prctile(Mg_fit1,prc_intervals(i),2); flipud(prctile(Mg_fit1,prc_intervals(end+1-i),2))],cmap(i,:),'edgecolor','none');
end
for i = 1:size(Antonelli_Mg,1)
    plot([Antonelli_Mg(i,1),Antonelli_Mg(i,1)],[Antonelli_Mg(i,2),Antonelli_Mg(i,4)],'-','linewidth',4,'color',[1 1 1]/3);
end
pM = plot(0,52,'kp','markersize',14,'markerfacecolor',[1 1 1]/3);
% for i = 1:num_iter
%     plot(xpred,Mg_fit1(:,i),'k-');
% end

set(gca,'xdir','reverse');
ylabel('[Mg] [mmol/kg]');
xlim(xlims);
ylim([0 60]);

%% Save Mg and Ca prediction curves

prc_intervals = 10:10:90;
Mg_pred = nan(length(xpred),length(prc_intervals));
Ca_pred = nan(length(xpred),length(prc_intervals));
for i = 1:length(prc_intervals)
    Mg_pred(:,i) = prctile(Mg_fit1,prc_intervals(i),2);
    Ca_pred(:,i) = prctile(Ca_fit1,prc_intervals(i),2);
end

save('MgCa_fit.mat','prc_intervals','xpred','Mg_pred','Ca_pred','Ca_fit1','Mg_fit1');
