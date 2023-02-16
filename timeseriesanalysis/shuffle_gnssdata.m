function [tobsfit,Dobsfit,Dobsfitn,tobspred,Dobspred,Dobspredn,boxYmin,boxYmax,Rmin,Rmax,dstd,dstdn]=shuffle_gnssdata(tobs,Dobs,tsfit,tspred)
% Author: Sambuddha Dhar

tobsfit=[];Dobsfit=[]; Rmin=[]; Rmax=[];
for i=1:size(Dobs,1)
    id=(tobs(i,:)>= min(tsfit) & tobs(i,:)<= max(tsfit));
    tobsfit(i,:)=tobs(i,id);
    Dobsfit(i,:)=Dobs(i,id);
    [Rmin(i,:),Rmax(i,:)]=deal(0,1);
    boxYmin(i,:)=min(Dobs(i,id));
    boxYmax(i,:)=max(Dobs(i,id));
end

Dobsfitn=Rmin+((Dobsfit-boxYmin).*(Rmax-Rmin)./(boxYmax-boxYmin));

for i=1:size(Dobs,1)
    id=(tobs(i,:)>= min(tspred) & tobs(i,:)<= max(tspred));
    tobspred(i,:)=tobs(i,id);
    Dobspred(i,:)=Dobs(i,id);
end

Dobspredn=(Dobspred-boxYmin)./(boxYmax-boxYmin);

stndev=@(x)  sqrt( sum((x-movmean(x,180,2)).^2, 2,'omitnan')./(size(x,2)-1) );
dstd=stndev(Dobspred);
dstdn=dstd./(boxYmax-boxYmin);
end