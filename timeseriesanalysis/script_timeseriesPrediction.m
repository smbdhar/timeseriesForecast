%% Main script for time series prediction
% Author: Sambuddha Dhar

% Related article: 
%   Function Model Based on Nonlinear Transient Rheology of Rocks:
%   An Analysis of Decadal GNSS Time Series After the 2011 Tohoku-oki Earthquake
%   Sambuddha Dhar and Jun Muto

% Pick any GNSS station's data set and run the model to forecast

clear;clc;
gnssid='960549'; stname={'Yamoto'}; % [T] % Put the data set here

gg1=importdata(['./',gnssid,'.mat']); stid={gnssid};
gg1stname={[stname{1},' NS'],[stname{1},' EW'],[stname{1},' Up']}';         
% ------------------------------------------------------------------- %
raw=[{gg1}];
[tobs,Dobs,Dint]=assemble_gnssdata(raw); %tobs=tobs(1,:); Dobs=Dobs(1,:); 
% tsfit=linspace(0.00,2,731);
tsfit=linspace(0,3.9,1424);
tspred=linspace(0,10,3651);
% tspred=linspace(0,5,1826);
[tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd]=shuffle_gnssdata(tobs,Dobs,tsfit,tspred);
Ndobs=size(Dobsfit,1);

% results (recent)
[theta1,theta2,theta3,theta4]=deal(0.0427527108314662,0.975363685964396,0.00300442746400926,0.968378142612612);

% 
[residual,Dcalpred,VRcalpred,AScalpred,cval]=linearRegression_gps_fixtheta(theta1,theta2,theta3,theta4,tsfit,tspred,...
    tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd);
obsft=Dobsfit;obspd=Dobspred;VRASpd=Dcalpred; VRpd=VRcalpred;ASpd=AScalpred; 


figure(9999);clf;set(gcf,'name','TS pickup and match');
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

obsflt=@(Dinp) [movmean(Dinp(1:7),1),movmean(Dinp(8:12),7),movmean(Dinp(13:end),14)];
obsfltw=@(Dinp) [movmean(Dinp(1:7),14),movmean(Dinp(8:12),14),movmean(Dinp(13:end),14)];

subplot(1,3,1); hold on; 
        plot(Ntobspred,Nobspd,'k*'); % Obs. pred
        p2=plot(tobsfit(1,:),obsft(1,:)+vssfitN,'g*'); p2.Color = [255 215 0]/255;% Obs. fit
        plot(Ntobspred,obsflt(Nobspd),'g'); %
        plot(Ntobspred,vssN,'g--');
        plot(tspred,NVRASpd,'r','LineWidth',2);        %
        plot(tspred,NVRpd,'b',tspred,NASpd,'g',tspred,vssVRASN,'g--'); 
        title(sitenm{1},'FontSize',18)%%
        xlabel('Years since earthquake')
        ylabel('disp. (m)'); box on;grid on;
subplot(1,3,2); hold on; 
        plot(Etobspred,Eobspd,'k*');            %
        p2=plot(tobsfit(2,:),obsft(2,:)+vssfitE,'g*'); p2.Color = [255 215 0]/255; %
        plot(Etobspred,obsflt(Eobspd),'g'); %
        plot(Etobspred,vssE,'g--');
        plot(tspred,EVRASpd,'r','LineWidth',2);        %
        plot(tspred,EVRpd,'b',tspred,EASpd,'g');                        %
        title(sitenm {2},'FontSize',18)%%
        xlabel('Years since earthquake')
        ylabel('disp. (m)'); box on;grid on;
subplot(1,3,3); hold on; 
        plot(Utobspred,Uobspd,'k*');            %
        p2=plot(tobsfit(3,:),obsft(3,:)+vssfitU,'g*'); p2.Color = [255 215 0]/255; %
        plot(Utobspred,obsfltw(Uobspd),'g'); %
        plot(Utobspred,vssU,'g--');
        plot(tspred,UVRASpd,'r','LineWidth',2);        %
        plot(tspred,UVRpd,'b',tspred,UASpd,'g'); 
        title(sitenm {3},'FontSize',18)%%
        xlabel('Years since earthquake')
        ylabel('disp. (m)'); box on;grid on;
        
% Error Envelop generation 
div=5;   
[X,Y,Z]=meshgrid([linspace(0.0182,0.1412,div),theta1],...
    [linspace(0.2846,3.2883,div),theta2],...
    [linspace(0.0026,0.0040,div),theta3]); 

XYZ=[X(:),Y(:),Z(:)];        
VRASsuit=repmat({zeros(length(tspred),size(XYZ,1))},1,Ndobs); VRsuit=VRASsuit; ASsuit=VRASsuit;
for ii=1:size(XYZ,1)
    [theta1,theta2,theta3]=deal(XYZ(ii,1),XYZ(ii,2),XYZ(ii,3));
    [~,Dcalpred,VRcalpred,AScalpred]=linearRegression_gps_fixtheta(theta1,theta2,theta3,theta4,tsfit,tspred,...
    tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd);
    k=[1:Ndobs];  
    obsft=Dobsfit;obspd=Dobspred;VRASpd=Dcalpred; VRpd=VRcalpred; ASpd=AScalpred;
    for i=1:length(k)  
        VRASsuit{i}(:,ii)=VRASpd(k(i),:)';%
        VRsuit{i}(:,ii)=VRpd(k(i),:)';%
        ASsuit{i}(:,ii)=ASpd(k(i),:)';%
    end
    clear('obspd','obsft','VRASpd','VRpd','ASpd'); 
end
VRASsuit{1}=VRASsuit{1}+vssVRASN(:);
VRASsuit{2}=VRASsuit{2}+vssVRASE(:);
VRASsuit{3}=VRASsuit{3}+vssVRASU(:);
if 1  %<***************************************************** %
    k=[1:Ndobs];  minErr=zeros(3,length(tspred)); maxErr=minErr;
    figure(9999);hold on;     
    obsft=Dobsfit;obspd=Dobspred;VRASpd=Dcalpred; envX=[];envY=[]; 
    for i=1:length(k)
        subplot(1,3,i); hold on;
        VRASenvmin=min(VRASsuit{i},[],2); VRASenvmax=max(VRASsuit{i},[],2);
        minErr(i,:)=VRASenvmin; maxErr(i,:)=VRASenvmax;
        mtx=[tspred,fliplr(tspred)]; mty=[VRASenvmin',fliplr(VRASenvmax')];
        patch(mtx,mty,'r','EdgeColor','none')
        envX(:,i)=mtx; envY(:,i)=mty;     
    end
    clear('obspd','obsft','VRASpd','VRpd','ASpd'); 
end

%% Residual plot
Nres=Nobspd-interp1(tspred,NVRASpd,Ntobspred);
Nminres=Nobspd-interp1(tspred,envY([1:length(tspred)],1)',Ntobspred);
Nmaxres=Nobspd-interp1(tspred,flipud(envY([length(tspred)+1:end],1))',Ntobspred);
NenvX=[Ntobspred,fliplr(Ntobspred)]; NenvY=[Nminres,fliplr(Nmaxres)];

Eres=Eobspd-interp1(tspred,EVRASpd,Etobspred);
Eminres=Eobspd-interp1(tspred,envY([1:length(tspred)],2)',Etobspred);
Emaxres=Eobspd-interp1(tspred,flipud(envY([length(tspred)+1:end],2))',Etobspred);
EenvX=[Etobspred,fliplr(Etobspred)]; EenvY=[Eminres,fliplr(Emaxres)];

Ures=Uobspd-interp1(tspred,UVRASpd,Utobspred);
Uminres=Uobspd-interp1(tspred,envY([1:length(tspred)],3)',Utobspred);
Umaxres=Uobspd-interp1(tspred,flipud(envY([length(tspred)+1:end],3))',Utobspred);
UenvX=[Utobspred,fliplr(Utobspred)]; UenvY=[Uminres,fliplr(Umaxres)];

figure(99992);clf;
subplot(1,3,1);hold on;
    patch(NenvX,NenvY,'g','EdgeColor','none');
    plot(Ntobspred,Nres,'r');
    plot(Ntobspred,Nminres,'b');
    plot(Ntobspred,Nmaxres,'c');
subplot(1,3,2);hold on;
    patch(EenvX,EenvY,'g','EdgeColor','none');
    plot(Etobspred,Eres,'r');
    plot(Etobspred,Eminres,'b');
    plot(Etobspred,Emaxres,'c');
subplot(1,3,3);hold on;
    patch(UenvX,UenvY,'g','EdgeColor','none');
    plot(Utobspred,Ures,'r');
    plot(Utobspred,Uminres,'b');
    plot(Utobspred,Umaxres,'c');    
% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = %