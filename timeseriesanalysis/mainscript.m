%% Main script for function model
% Author: Sambuddha Dhar

% Related article: 
%   Function Model Based on Nonlinear Transient Rheology of Rocks:
%   An Analysis of Decadal GNSS Time Series After the 2011 Tohoku-oki Earthquake
%   Sambuddha Dhar and Jun Muto

% This script is used for training the function model using six GNSS data
% set; Same is also used to predict time series for other stations.
% The GNSS time series are publicly available at http://terras.gsi.go.jp/
% upon free registrations. 
% For each dataset, .mat file = [t,N,E,U, VstN, VstE, VstU]
% where t=time; (N)orth (E)ast (U)p and Vst stands for steady-state disp. 

clear;clc;
gg1=importdata('./dataset/960549.mat'); % Yamoto
gg2=importdata('./dataset/940028.mat'); % Miyako
gg3=importdata('./dataset/940042.mat'); % Iwaki
gg4=importdata('./dataset/950193.mat');% Minase
sitename={'Yamoto 960549 NS','Yamoto 960549 EW','Yamoto 960549 Up',...
          'Miyako 940028 NS','Miyako 940028 EW','Miyako 940028 Up',...
          'Iwaki 940041 NS','Iwaki 940041 EW','Iwaki 940041 Up ',...
          'Minase 950193 NS','Minase 950193 EW','Minase 950193 Up'};
raw=[{gg1};{gg2};{gg3};{gg4}];%
gg5=importdata('./dataset/950233.mat'); % 950232 <> 960565
    gg5stname={'Ryoutsu2 950233 NS','Ryoutsu2 950233 EW','Ryoutsu2 950233 Up'};
gg6=importdata('./dataset/950151.mat'); % Kanita 950151 <>950182
    gg6stname={'Kanita 950151 NS','Kanita 950151 EW','Kanita 950151 Up'}; 

raw=[raw;gg5;gg6]; sitename=[sitename,gg5stname,gg6stname]; % ;gg7 ,gg7stname
siteID={'960549','940028','940042','950193','950233','950151'}; %,'940027'

[tobs,Dobs,Dint]=assemble_gnssdata(raw); %
% tsfit=linspace(0.0,1,366);
% tsfit=linspace(0.0,2,731);
tsfit=linspace(0,3.9,1424);
% tsfit=linspace(0.0,5.8,2118);
% tsfit=linspace(0.0,8.9,3250);
 
tspred=linspace(0,10,3651);
[tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd,dstdn]=shuffle_gnssdata(tobs,Dobs,tsfit,tspred);
Ndobs=size(Dobsfit,1);
%% Fitting and ploting observed and modeled time series
[theta1,theta2,theta3,theta4]=deal(0.0427527108314662,0.975363685964396,0.00300442746400926,0.968378142612612);

[residual,Dcalpred,VRcalpred,AScalpred]=linearRegression_gnss_theta(theta1,theta2,theta3,theta4,tsfit,tspred,...
    tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd);
residual

% ************* plot timeseries *****************************%%
% if 1  %<***************************************************** %                    %
    figure(991);clf;set(gcf,'name','TS obs+fit view')          %                                                     %
        obsfit=Dobsfit;obspd=Dobspred;VRASpd=Dcalpred;         %
        VRpd=VRcalpred;ASpd=AScalpred; 
        nots=size(obsfit,1)./3;
        last=4;
%         kid=[1,2,3,4];
for     k=1:last
  for   kk=[3.*(k-1)+1, 3.*(k-1)+2, 3.*(k-1)+3]
      
        vss=interp1(tobs(kk,:),Dint(kk,:),tobspred(kk,:));
        vssfit=interp1(tobs(kk,:),Dint(kk,:),tobsfit(kk,:));
        vssVRAS=interp1(tobs(kk,:),Dint(kk,:),tspred);
        
        subplot(last,3,kk); hold on; 
        plot(tobspred(kk,:),obspd(kk,:)+vss,'k*');
        plot(tobspred(kk,:),vss,'g--');
        p2=plot(tobsfit(kk,:),obsfit(kk,:)+vssfit,'k*');p2.Color = [255 215 0]/255; 
        plot(tspred,VRpd(kk,:),'b',tspred,ASpd(kk,:),'g','LineWidth',2);
        plot(tspred,VRASpd(kk,:)+vssVRAS,'r','LineWidth',2); 
        title(sitename{kk},'FontSize',18)%%
        xlabel('Years since earthquake')
        ylabel('disp. (m)'); box on;grid on;
  end     
end

if 0  % ***** plot timeseries residual************************ %
    figure(992);clf;set(gcf,'name','TS Residual')          %                                                     %
    for     k=1:last
        for kk=[3.*(k-1)+1, 3.*(k-1)+2, 3.*(k-1)+3]
            subplot(last,3,kk); hold on;
            res=obspd(kk,:)-interp1(tspred,VRASpd(kk,:),tobspred(kk,:),'linear');
            plot(tobspred(kk,:),res,'k*');
            title(sitename{kk},'FontSize',18)%%
            ylim([-0.06 0.06])
            xlabel('Years since earthquake');ylabel('disp. (m)'); box on;grid on;
        end
    end
end
% ************************************************************ %
%% For Bayesian search using GP-UCB method
d=4; % Dimension
N=1000; % Initial sampling (nodal points)
n=10000; % Test points number (not related to sample)
iter=900;
bounds=[0.001 10;...% 
        0.001 10;...% 
        0.00033 0.0066;... % lambda=actual lambda/300 MPa
        0.001 1000]; %4.5467 0.001 1000
f=@(m) linearRegression_gnss_theta(m(1),m(2),m(3),m(4),tsfit,tspred,...
    tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd);
%
tic;
[bestX,bestscore,X,y]=bayesOptlog(f,bounds,N,d,n,iter,'plot_true'); %
toc;
bestX
% save('bestX.mat','bestX');