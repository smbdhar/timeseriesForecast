function [tobs,Dobs,Dint]=assemble_gnssdata(raw)
% Author: Sambuddha Dhar

% raw=[{gg1};{gg2};{gg3};{gg4}]; %cell data
delT= year2decimal(2011,03,12,09,00,01)-year2decimal(2011,03,11,14,56,00);

tt=[];
for i=1:length(raw)
    tt=[tt,raw{i}(:,1)'];
end
tt=unique(tt);
tobs=[]; Dobs=[]; Dint=[];
for i=1:length(raw)
    tdata=tt(:)';
    north=interp1(raw{i}(:,1),raw{i}(:,2),tt,'linear');
    east=interp1(raw{i}(:,1),raw{i}(:,3),tt,'linear');
    up=interp1(raw{i}(:,1),raw{i}(:,4),tt,'linear');

    Intnorth=interp1(raw{i}(:,1),raw{i}(:,5),tt,'linear');
    Inteast=interp1(raw{i}(:,1),raw{i}(:,6),tt,'linear');
    Intup=interp1(raw{i}(:,1),raw{i}(:,7),tt,'linear');
    
    id=find(tdata>=delT);
    tdata=tdata(id)-delT;
    north=north(id);
    east=east(id);
    up=up(id);
    Intnorth=Intnorth(id);
    Inteast=Inteast(id);
    Intup=Intup(id);
    
    north=north-north(1);
    east=east-east(1);
    up=up-up(1);
    Intnorth=Intnorth-Intnorth(1);
    Inteast=Inteast-Inteast(1);
    Intup=Intup-Intup(1);
    
    tobs=[tobs;tdata;tdata;tdata]; Dobs=[Dobs;north(:)';east(:)';up(:)'];
    Dint=[Dint;Intnorth(:)';Inteast(:)';Intup(:)'];
end
% tobs=tobs-tt(1);
end