function [L2err,lmcosi]=fitOoid(v,L)
% [L2err,lmcosi]=ooid(v,L)
%
% EXAMPLE:
%
% for index=1:10
%   L2err(index)=ooid(v,index); pause
% end
    
%load v
x=v(:,1);
y=v(:,2);
z=v(:,3);

% Take out the center
x=x-mean(x);
y=y-mean(y);
z=z-mean(z);

% My definitions
[phi,th,robs]=cart2sph(x,y,z);
% Latitude in degrees
lat=th/pi*180;
% Longitude in degrees
lon=phi*180/pi;

% Truncation
%L=2; % Ellipsoidal

% Do not worry about the first three warnings

% Expansion 
[lmcosi,~,L2err]=xyz2plm(robs,L,'irr',lat,lon);
%%
figure(1); clf
%tiledFigure(2.5, 2.5)
%subplot(121)
s = scatter3(x,y,z,10,'filled', 'k');
s.MarkerFaceAlpha = .5;
s.MarkerEdgeColor = 'none';

axis equal on
ax=axis;
rh=plm2xyz(lmcosi,lat,lon);
hold on
[xh,yh,zh]=sph2cart(lon*pi/180,lat*pi/180,rh);
scatter3(xh,yh,zh,20, [.8 0 0])
set(gca,'XTickLabels', [],'YTickLabels', [],'ZTickLabels', [])
hold off
grid on
box on
%title(sprintf('approximation error %g',L2err))

%%
%keyboard
% subplot(223)
% plot([x ; y ; z],[xh ; yh ; zh],'+')
% axis equal
% grid on
%{
subplot(122)
[r,lon,lat]=plm2xyz(lmcosi,1);
% Map
%plotplm(lmcosi,[],[],[],1)
plotonsphere(r,0.2)
axis equal on
%axis(ax)
grid on

title(sprintf('L = %i approximation',L))
%}
end