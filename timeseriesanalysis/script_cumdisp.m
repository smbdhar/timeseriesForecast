%% Main script for calculating cumulative displacements
% Author: Sambuddha Dhar
% Related article: 
%   Function Model Based on Nonlinear Transient Rheology of Rocks:
%   An Analysis of Decadal GNSS Time Series After the 2011 Tohoku-oki Earthquake
%   Sambuddha Dhar and Jun Muto

% Cumulative displacements- Obs, model, VR, AS, residuals

clear;clc;
[stid,lon,lat,sitename]=textread('./networkNEJP.dat','%s %s %s %s','commentstyle','shell');
% [stid,lon,lat,sitename]=textread('./networkNEJP_train.dat','%s %s %s %s','commentstyle','shell');

path1='./dataset';
path2='./dataset';

timept=[0.5,1,2,3,4,5,6,7,8,9,10-(2/365)]; %for residuals

Nobsmat=zeros(length(stid),length(timept)); NVssmat=zeros(length(stid),length(timept)); 
Nresmat=Nobsmat; NVRASmat=Nobsmat; NVRmat=Nobsmat; NASmat=Nobsmat;  
Eobsmat=zeros(length(stid),length(timept)); EVssmat=zeros(length(stid),length(timept));
Eresmat=Eobsmat; EVRASmat=Eobsmat; EVRmat=Eobsmat; EASmat=Eobsmat;
Uobsmat=zeros(length(stid),length(timept)); UVssmat=zeros(length(stid),length(timept));
Uresmat=Uobsmat; UVRASmat=Uobsmat; UVRmat=Uobsmat; UASmat=Uobsmat;
%%
for i=1:length(stid)
    gnssid=stid{i}; stname=sitename(i);
    
filepath1=[path1,'/',gnssid,'.mat'];
filepath2=[path2,'/',gnssid,'.mat'];
if exist(filepath1)
    filepath=filepath1;
else 
    filepath=filepath2;
end
gg1=importdata(filepath); %stid={gnssid};
gg1stname={[stname{1},' NS'],[stname{1},' EW'],[stname{1},' Up']}';         
%% ------------------------------------------------------------------- %
raw=[{gg1}];
% sitename=[gg1stname];
[tobs,Dobs,Dint]=assemble_gnssdata(raw); %tobs=tobs(1,:); Dobs=Dobs(1,:);
% tsfit=linspace(0.0,1,366);
% tsfit=linspace(0.00,2,731);
% tsfit=linspace(0,3.9,1424);
% tsfit=linspace(0.0,5.8,2118);
 tsfit=linspace(0.0,8.9,3250);

tspred=linspace(0,10,3651);
[tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd]=shuffle_gnssdata(tobs,Dobs,tsfit,tspred);
Ndobs=size(Dobsfit,1);

%t=1
% [theta1,theta2,theta3,theta4]=deal(0.0385833224643496,0.532728166620938,0.00205340378058017,30.3680625880763);
%t=2
% [theta1,theta2,theta3,theta4]=deal(0.0860951489476105,0.741724447809488,0.00271351338129483,1.29594523852338);
%t=5.8
% [theta1,theta2,theta3,theta4]=deal(0.0869689187148345,1.47627764637125,0.00310037846255843,0.313663007224935);
%t=8.9
% [theta1,theta2,theta3,theta4]=deal(1.23345552686723,0.0794013415888254,0.00309920021027875,0.244504055832404);


% best model
[theta1,theta2,theta3,theta4]=deal(0.0427527108314662,0.975363685964396,0.00300442746400926,0.968378142612612);


% ==== ==== ==== ===== 
[residual,Dcalpred,VRcalpred,AScalpred,cval]=linearRegression_gnss_theta(theta1,theta2,theta3,theta4,tsfit,tspred,...
    tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd);
obsft=Dobsfit;obspd=Dobspred;VRASpd=Dcalpred; VRpd=VRcalpred;ASpd=AScalpred; 

% figure(9999);clf;set(gcf,'name','TS pickup and match');
    VN=(Dint(1,2)-Dint(1,1))./(tobs(1,2)-tobs(1,1));
    VE=(Dint(2,2)-Dint(2,1))./(tobs(2,2)-tobs(2,1));
    VU=(Dint(3,2)-Dint(3,1))./(tobs(3,2)-tobs(3,1));
     
    vssN=tobspred(1,:).*VN; vssfitN=tobsfit(1,:).*VN; vssVRASN=tspred.*VN;               
    vssE=tobspred(2,:).*VE; vssfitE=tobsfit(2,:).*VE; vssVRASE=tspred.*VE;       
    vssU=tobspred(3,:).*VU; vssfitU=tobsfit(3,:).*VU; vssVRASU=tspred.*VU;
        
Ntobspred=tobspred(1,:); Etobspred=tobspred(2,:); Utobspred=tobspred(3,:);
Nobspd=obspd(1,:)+vssN; Eobspd=obspd(2,:)+vssE; Uobspd=obspd(3,:)+vssU;
NVRASpd=VRASpd(1,:)+vssVRASN; EVRASpd=VRASpd(2,:)+vssVRASE; UVRASpd=VRASpd(3,:)+vssVRASU;
NVRpd=VRpd(1,:); EVRpd=VRpd(2,:); UVRpd=VRpd(3,:);
NASpd=ASpd(1,:); EASpd=ASpd(2,:); UASpd=ASpd(3,:);
sitenm=gg1stname;%[sitename{1}];

% calculating cumulative displacemements
obsflt=@(Dinp) [movmean(Dinp(1:7),1),movmean(Dinp(8:12),14),movmean(Dinp(13:end),14)];
obsfltw=@(Dinp) [movmean(Dinp(1:7),14),movmean(Dinp(8:12),14),movmean(Dinp(13:end),14)];

NcumObs=obsflt(Nobspd); NcumObs=NcumObs-NcumObs(1);
NcumVRAS=interp1(tspred,NVRASpd,Ntobspred); %NcumVRAS=NcumVRAS-NcumVRAS(1);
NcumVR=interp1(tspred,NVRpd,Ntobspred); %NcumVR=NcumVR-NcumVR(1);
NcumAS=interp1(tspred,NASpd,Ntobspred); %NcumAS=NcumAS-NcumAS(1);

EcumObs=obsflt(Eobspd); EcumObs=EcumObs-EcumObs(1);
EcumVRAS=interp1(tspred,EVRASpd,Etobspred); %EcumVRAS=EcumVRAS-EcumVRAS(1);
EcumVR=interp1(tspred,EVRpd,Etobspred); %EcumVR=EcumVR-EcumVR(1);
EcumAS=interp1(tspred,EASpd,Etobspred); %EcumAS=EcumAS-EcumAS(1);

UcumObs=obsfltw(Uobspd); UcumObs=UcumObs-UcumObs(1);
UcumVRAS=interp1(tspred,UVRASpd,Utobspred); %UcumVRAS=UcumVRAS-UcumVRAS(1);
UcumVR=interp1(tspred,UVRpd,Utobspred); %UcumVR=UcumVR-UcumVR(1);
UcumAS=interp1(tspred,UASpd,Utobspred); %UcumAS=UcumAS-UcumAS(1);

% snapshot of cumulative displacements
NcumObspt=interp1(Ntobspred,NcumObs,[timept]);  %(14/365) 
%     if isnan(NcumObspt(1)); NcumObspt(1)=0;end
%     NcumObspt=NcumObspt(2:end)-NcumObspt(1);
NcumVRASpt=interp1(Ntobspred,NcumVRAS,[timept]); %(14/365) 
%     if isnan(NcumVRASpt(1)); NcumVRASpt(1)=0;end
%     NcumVRASpt=NcumVRASpt(2:end)-NcumVRASpt(1);

NcumVRpt=interp1(Ntobspred,NcumVR,timept);
NcumASpt=interp1(Ntobspred,NcumAS,timept);  
    Ncumvsspt=interp1(Ntobspred,vssN,timept);

EcumObspt=interp1(Etobspred,EcumObs,[timept]); %(14/365) 
%     if isnan(EcumObspt(1)); EcumObspt(1)=0;end
%     EcumObspt=EcumObspt(2:end)-EcumObspt(1);
EcumVRASpt=interp1(Etobspred,EcumVRAS,[timept]); %(14/365) 
%     if isnan(EcumVRASpt(1)); EcumVRASpt(1)=0;end
%     EcumVRASpt=EcumVRASpt(2:end)-EcumVRASpt(1);
EcumVRpt=interp1(Etobspred,EcumVR,timept);
EcumASpt=interp1(Etobspred,EcumAS,timept);
    Ecumvsspt=interp1(Etobspred,vssE,timept);
    
UcumObspt=interp1(Utobspred,UcumObs,[ timept]);
%     if isnan(UcumObspt(1)); UcumObspt(1)=0;end
%     UcumObspt=UcumObspt(2:end)-UcumObspt(1);
UcumVRASpt=interp1(Utobspred,UcumVRAS,[ timept]); 
%     if isnan(UcumObspt(1)); UcumObspt(1)=0;end
%     UcumVRASpt=UcumVRASpt(2:end)-UcumVRASpt(1);
UcumVRpt=interp1(Utobspred,UcumVR,timept);
UcumASpt=interp1(Utobspred,UcumAS,timept);
    Ucumvsspt=interp1(Utobspred,vssU,timept);

Nobsmat(i,:)=NcumObspt; Nresmat(i,:)=NcumObspt-NcumVRASpt; NVssmat(i,:)=Ncumvsspt;%VN;
NVRASmat(i,:)=NcumVRASpt;
NVRmat(i,:)=NcumVRpt;
NASmat(i,:)=NcumASpt;

Eobsmat(i,:)=EcumObspt; Eresmat(i,:)=EcumObspt-EcumVRASpt; EVssmat(i,:)=Ecumvsspt;%VE;
EVRASmat(i,:)=EcumVRASpt;
EVRmat(i,:)=EcumVRpt;
EASmat(i,:)=EcumASpt;

Uobsmat(i,:)=UcumObspt; Uresmat(i,:)=UcumObspt-UcumVRASpt; UVssmat(i,:)=Ucumvsspt;%VU;
UVRASmat(i,:)=UcumVRASpt;
UVRmat(i,:)=UcumVRpt;
UASmat(i,:)=UcumASpt;
end
