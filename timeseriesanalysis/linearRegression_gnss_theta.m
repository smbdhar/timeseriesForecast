function [varargout]=linearRegression_gnss_theta(theta1,theta2,theta3,theta4,tsfit,tspred,tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd) 
% Author: Sambuddha Dhar
% 
Ndobs=size(Dobsfit,1);
Dcalpred=zeros(Ndobs,length(tspred));
VRcalpred=zeros(Ndobs,length(tspred));
AScalpred=zeros(Ndobs,length(tspred));
resid=zeros(Ndobs,1); cval=zeros(Ndobs,3);
for kk=1:Ndobs

tobsfitk=tobsfit(kk,:);Dobsfitk=Dobsfit(kk,:);Dobsfitnk=Dobsfitn(kk,:);
tobspredk=tobspred(kk,:);Dobspredk=Dobspred(kk,:);Dobsprednk=Dobspredn(kk,:);
boxYmink=boxYmin(kk,:);boxYmaxk=boxYmax(kk,:);%[Rmink,Rmaxk]=deal(Rmin(kk,:),Rmax(kk,:));
dstdk=dstd(kk,:);

[error,mb]=linearReg(theta1,theta2,theta3,tobsfitk,Dobsfitk,dstdk,tsfit);
cval(kk,:)=mb;
% figure(9);clf;hold on;
% plot(tobspredk,Dobspredk,'k*');plot(tobsfitk,Dobsfitk,'g*');
% [dcal,VR,AS,bestm]=anlymodel(tsfit,[mb,theta1,theta2,theta3],'none');
% [Dcalprednk,VRcalprednk,AScalprednk]=anlymodel(tspred,[bestm,theta1,theta2,theta3],'none');
[Dcalpredk,VRcalpredk,AScalpredk]=anlyModel(tspred,[mb,theta1,theta2,theta3],'none');

% Dcalpredk=boxYmink+((Dcalprednk-Rmink).*(boxYmaxk-boxYmink)./(Rmaxk-Rmink));
% VRcalpredk=boxYmink+((VRcalprednk-Rmink).*(boxYmaxk-boxYmink)./(Rmaxk-Rmink));
% AScalpredk=boxYmink+((AScalprednk-Rmink).*(boxYmaxk-boxYmink)./(Rmaxk-Rmink));
    Dcalpredm=interp1(tspred,Dcalpredk,tobspredk);
    
   obsflt=@(Dinp) [movmean(Dinp(1:7),2),movmean(Dinp(8:12),14),movmean(Dinp(13:end),14)];
%     obsfltw=@(Dinp) [movmean(Dinp(1:7),14),movmean(Dinp(8:12),14),movmean(Dinp(13:end),14)];

   Dobspredk=obsflt(Dobspredk);
    
    res=sum((Dobspredk-Dcalpredm).^2,2,'omitnan');%./(dstdk.^2);
Dcalpred(kk,:)=Dcalpredk;
VRcalpred(kk,:)=VRcalpredk;
AScalpred(kk,:)=AScalpredk;
resid(kk,:)=res;
end
residual=-sum(resid,1);%./numel(Dobspred);
residual=1e3.*exp(0.5.*residual./theta4);
residual=((2*pi*theta4).^(-3/2)).*residual;
% plot(tspred,Dcalpredk,'r')
varargout{1}=residual; 
varargout{2}=Dcalpred; 
varargout{3}=VRcalpred; 
varargout{4}=AScalpred;
varargout{5}=cval;
end

function [error,mb]=linearReg(theta1,theta2,theta3,tobs,dobs,dstd,ts)
f=@(m) lossfunc([m(1),m(2),m(3),theta1,theta2,theta3],tobs,dobs,ts);
lb=[-Inf,-Inf,-Inf];%0-dstd
ub=[+Inf,+Inf,+Inf]; %1+dstd
m0=rand(1,3); %(lb+ub)/2;%
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500,'MaxFunctionEvaluations',2000,'Display','off');
mb = fmincon(f,m0,[],[],[],[],lb,ub,[],options);

% options=optimset('Display','off','TolX',1e-16, 'MaxIter',1500,'MaxFunEvals',2000);
% mb = fminsearch(f,m0,options);
% opts = optimset('fminsearch');
% opts.TolFun = 1.e-12;

% opts = optimset('Display','off','MaxFunEvals',500000, 'MaxIter',100000,'TolFun',1e-30);
% [mb] = fminsearchbnd(f,m0,lb,ub,opts);

error=f(mb);
function [error]=lossfunc(m,tobs,dobs,ts)
% m=row; t=row; tobs=[1 x tobs]; dobs=[1 x tobs];
ts=ts(:)';
% dobsmu=movmean(dobs,5);
dobsmu=dobs;

obsC=interp1(tobs,dobsmu,ts);
[dc]=anlyModel(ts,m,'none');
res=(obsC-dc).^2;
error=sum(res,2,'omitnan');

% dc1=interp1(ts,dc,tobs); deldc1=dc1(end)-dc1(1);
% deldobsmu=dobsmu(end)-dobsmu(1);
% error=(deldobsmu-deldc1).^2;

end
end
%