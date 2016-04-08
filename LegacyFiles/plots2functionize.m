VarList={'LL04','RL04','LL10','RL10'};

varIDX =4;
PC = 1;
figure;
[sorted,idx] = sort(data.(VarList{varIDX}).PCA.score(:,PC));
imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);

%% this is the basic 4 way PCA plot. this could be adapted to be useful for lots of things.
VarList={'LL04','RL04','LL10','RL10'};
PC = 1;
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    [sorted,idx] = sort(data.(VarList{varIDX}).PCA.score(:,PC));
    imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);
    title(['PC: ' num2str(PC) '   ' VarList{varIDX}]);
end
%% get the zscore data for outputting for stats
for varIDX=1:4
    zDataOutput{varIDX}=zscore(data.(VarList{varIDX}).MlMtx')';
end
%% this is the overall means, zscored - plotted as 4 subplots
clearvars ('-except', 'data');
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;

for varIDX=1:4
    subplot(2,2,varIDX);
    zData=zscore(data.(VarList{varIDX}).MlMtx')';
    y=nanmean(zData);
    yStdErr=nanstd(zData)/sqrt(size(zData,1));
    plot(data.xA,y,'linewidth',.8);
    hold on
    %plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','.','markersize',3,'linestyle','none');
    plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    plot(downsample(data.xA,1),downsample(y+yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    %plot(data.xA,y+yStdErr,'linestyle',':');
    title(['Overall zscore: ' VarList{varIDX}]);
    axis([-3 3 -1 1])
end
%% overal means zscored - plotted as 1 plot
clearvars ('-except', 'data');
figure;
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:4
    zData=zscore(data.(VarList{varIDX}).MlMtx')';
    y(:,varIDX)=nanmean(zData);
    yStdErr(:,varIDX)=nanstd(zData)/sqrt(size(zData,1));
end

    %subplot(2,1,condIDX);
    plot(data.xA,y,'linewidth',2.5);
    hold on
    %plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','.','markersize',3,'linestyle','none');
    plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    plot(downsample(data.xA,1),downsample(y+yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    %plot(data.xA,y+yStdErr,'linestyle',':');
    title(['Overall zscore']);
    axis([-3 3 -1 1]);
    legend (VarList{:});

%% custom sorting  0 -> .5s
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    [sorted, idx] = sort(mean(tempZ(:,data.xA>=0 & data.xA<=.5),2 ));
    imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);
    title(['0-.5 zscore sorted: ' VarList{varIDX}]);
end

%% custom sorting  0 -> 2.5s
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    [sorted, idx] = sort(mean(tempZ(:,data.xA>=0 & data.xA<=2.5),2 ));
    imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);
    title(['0-2.5 zscore sorted: ' VarList{varIDX}]);
end
%% custom sorting  0 -> 5s
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    [sorted, idx] = sort(mean(tempZ(:,data.xA>=0 & data.xA<=5),2 ));
    imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);
    title(['0-5 zscore sorted: ' VarList{varIDX}]);
    data.(VarList{varIDX}).mySort.idx=idx;
end

%% custom sorting  0 -> 5s -- reversing sorts
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:2
    subplot(2,2,varIDX*2);
    tempZ= zscore(data.(VarList{varIDX*2-1}).MlMtx')';
    [sorted, s1] = sort(mean(tempZ(:,data.xA>=0 & data.xA<=5),2 ));
    tempZ= zscore(data.(VarList{varIDX*2}).MlMtx')';
    [sorted, s2] = sort(mean(tempZ(:,data.xA>=0 & data.xA<=5),2 ));
    
    imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX*2-1}).MlMtx(s2,:)')'),[-2 2]);
    title(['0-5 zscore sorted: ' VarList{varIDX*2-1}]);
    subplot(2,2,varIDX*2-1);
    imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX*2}).MlMtx(s1,:)')'),[-2 2]);
    title(['0-5 zscore sorted: ' VarList{varIDX*2}]);
end
%% custom sorting 5.5-6.5
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    [sorted, idx] = sort(mean(tempZ(:,data.xA>=5.5 & data.xA<=6.5),2 ));
    imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);
    title(['5.5-6.5 zscore sorted: ' VarList{varIDX}]);
end

%% custom sorting diff -3to0 vs 0-3s
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    [sorted, idx] = sort(mean( tempZ(:,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(:,data.xA>=0 & data.xA<=3),2 ));
    imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);
    title(['-3>0 vs 0>3 zscore sorted: ' VarList{varIDX}]);
    data.(VarList{varIDX}).mySort.idx=idx;
end
%% Was also thinking that once you're good with the d prime analysis,  could compare latencies after press to change in sig modulated cells
%% Getting time between rf and bb

DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m22D04_03.mat';'m25D04_03.mat';...
    'm13D10_03.mat';'m14D04_03.mat';'m21D10_03.mat';'m22D10_03.mat';'m25D10_03.mat'};

clear tempDelay
for XX=1:length(DirList)
    for evtTrigger={'LL_C' 'RL_C'}
        load(DirList{XX});
        % Mean Firing Rates
        if strcmp(evtTrigger,'RL_C')
            Press_raw=strmatch('RL_R',behaveEvt_Raw);
            Press_ts=behaveEvtTm_Raw(Press_raw);
            BB_raw=strmatch('BB',behaveEvt_Raw);
            BB_ts=behaveEvtTm_Raw(BB_raw);
            tempDelay=arrayfun (@(x) BB_ts(find(BB_ts>x & BB_ts<x+4,1))-x, Press_ts, 'UniformOutput', false);
            %k=unique(cell2mat(arrayfun (@(x) find(BB_ts>x & BB_ts<x+8,1), Press_ts, 'UniformOutput', false)));
            delayTS(XX,2)=mean(cell2mat(tempDelay));
        elseif strcmp(evtTrigger,'LL_C')
            Press_raw=strmatch('LL_R',behaveEvt_Raw);
            Press_ts=behaveEvtTm_Raw(Press_raw);
            BB_raw=strmatch('BB',behaveEvt_Raw);
            BB_ts=behaveEvtTm_Raw(BB_raw);
            tempDelay=arrayfun (@(x) BB_ts(find(BB_ts>x & BB_ts<x+4,1))-x, Press_ts, 'UniformOutput', false);
            delayTS(XX,1)=mean(cell2mat(tempDelay));
        end
    end
end
%% Getting time until no more bb after rf to bb

DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m22D04_03.mat';'m25D04_03.mat';...
    'm13D10_03.mat';'m14D04_03.mat';'m21D10_03.mat';'m22D10_03.mat';'m25D10_03.mat'};

clear tempDelay
for XX=1:length(DirList)
    for evtTrigger={'LL_C' 'RL_C'}
        load(DirList{XX});
        % Mean Firing Rates
        if strcmp(evtTrigger,'RL_C')
            Press_raw=strmatch('RL_R',behaveEvt_Raw);
            Press_ts=behaveEvtTm_Raw(Press_raw);
            BB_raw=strmatch('BB',behaveEvt_Raw);
            BB_ts=behaveEvtTm_Raw(BB_raw);
            tempDelay=arrayfun (@(x) BB_ts(find(BB_ts>x & BB_ts<x+4,1))-x, Press_ts, 'UniformOutput', false);
            %k=unique(cell2mat(arrayfun (@(x) find(BB_ts>x & BB_ts<x+8,1), Press_ts, 'UniformOutput', false)));
            delayTS(XX,2)=mean(cell2mat(tempDelay));
        elseif strcmp(evtTrigger,'LL_C')
            Press_raw=strmatch('LL_R',behaveEvt_Raw);
            Press_ts=behaveEvtTm_Raw(Press_raw);
            BB_raw=strmatch('BB',behaveEvt_Raw);
            BB_ts=behaveEvtTm_Raw(BB_raw);
            tempDelay=arrayfun (@(x) BB_ts(find(BB_ts>x & BB_ts<x+4,1))-x, Press_ts, 'UniformOutput', false);
            delayTS(XX,1)=mean(cell2mat(tempDelay));
        end
    end
end
%% Set Inc and Decreasers from DPrime...This is for RFDelivery
VarList={'LL04','RL04','LL10','RL10'};
for varIDX=1:4
    tempZ=zscore(data.(VarList{varIDX}).MlMtx')';
    SIG=data.(VarList{varIDX}).DPrime.First3.SIG;
    CIDPIdx=data.(VarList{varIDX}).DPrime.First3.CIDPIdx;
    newIDs=sort(CIDPIdx(SIG));
    %     tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    %     [sorted, idx] = sort(mean( tempZ(:,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(:,data.xA>=0 & data.xA<=3),2 ));
    diffmean=mean( tempZ(newIDs,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(newIDs,data.xA>=0 & data.xA<=3),2 );
    data.(VarList{varIDX}).DPrime.First3.pieData(1,1)=length(find(diffmean<0));
    data.(VarList{varIDX}).DPrime.First3.pieData(1,3)=length(find(diffmean>0));
    data.(VarList{varIDX}).DPrime.First3.pieData(1,2)=size(tempZ,1)-size(SIG,1);
    data.(VarList{varIDX}).DPrime.First3.IncZ=tempZ(newIDs((diffmean<0)),:);
    data.(VarList{varIDX}).DPrime.First3.DecZ=tempZ(newIDs((diffmean>0)),:);
    data.(VarList{varIDX}).DPrime.First3.IncIDs=newIDs((diffmean<0));
    data.(VarList{varIDX}).DPrime.First3.DecIDs=newIDs((diffmean>0));
    data.(VarList{varIDX}).DPrime.First3.NCIDs=sort(CIDPIdx(1:SIG(1)-1));
    data.(VarList{varIDX}).DPrime.First3.NCZ=tempZ(data.(VarList{varIDX}).DPrime.First3.NCIDs);
    diffmean=mean( tempZ(CIDPIdx,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(CIDPIdx,data.xA>=0 & data.xA<=3),2 );
    diffmean=mean( tempZ(CIDPIdx,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(CIDPIdx,data.xA>=0 & data.xA<=3),2 );
    diffdecidx=find(diffmean<=0);
    diffincidx=find(diffmean>0);
    data.(VarList{varIDX}).DPrime.First3.overallSort=[flipud(diffincidx);diffdecidx];
end

%% Set Inc and Decreasers from DPrime...This is for RFCon
VarList={'LL04','RL04','LL10','RL10'};
for varIDX=1:4
    tempZ=zscore(data.(VarList{varIDX}).MlMtx')';
    SIG=data.(VarList{varIDX}).DPrime.First3.SIG;
    CIDPIdx=data.(VarList{varIDX}).DPrime.First3.CIDPIdx;
    newIDs=sort(CIDPIdx(SIG));
    %     tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    %     [sorted, idx] = sort(mean( tempZ(:,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(:,data.xA>=0 & data.xA<=3),2 ));
    diffmean=mean( tempZ(newIDs,data.xA>=0 & data.xA<=5),2)- mean(tempZ(newIDs,D),2 );
    data.(VarList{varIDX}).DPrime.First3.pieData(1,1)=length(find(diffmean<0));
    data.(VarList{varIDX}).DPrime.First3.pieData(1,3)=length(find(diffmean>0));
    data.(VarList{varIDX}).DPrime.First3.pieData(1,2)=size(tempZ,1)-size(SIG,1);
    data.(VarList{varIDX}).DPrime.First3.IncZ=tempZ(newIDs((diffmean<0)),:);
    data.(VarList{varIDX}).DPrime.First3.DecZ=tempZ(newIDs((diffmean>0)),:);
    data.(VarList{varIDX}).DPrime.First3.IncIDs=newIDs((diffmean<0));
    data.(VarList{varIDX}).DPrime.First3.DecIDs=newIDs((diffmean>0));
    data.(VarList{varIDX}).DPrime.First3.NCIDs=sort(CIDPIdx(1:SIG(1)-1));
    data.(VarList{varIDX}).DPrime.First3.NCZ=tempZ(data.(VarList{varIDX}).DPrime.First3.NCIDs);
    diffmean=mean( tempZ(CIDPIdx,data.xA>=0 & data.xA<=5),2)- mean(tempZ(CIDPIdx,data.xA>=8 & data.xA<=13),2 );
    diffdecidx=find(diffmean<=0);
    diffincidx=find(diffmean>0);
    data.(VarList{varIDX}).DPrime.First3.overallSort=[flipud(diffincidx);diffdecidx];
    figure; imagesc(flipud(tempZ(CIDPIdx(data.(VarList{varIDX}).DPrime.First3.overallSort),:)),[-2 2]);
end
%% DPrime Plot Inc and Dec
clearvars ('-except', 'data');
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:4
    zData=data.(VarList{varIDX}).DPrime.First3.IncZ;
    y(:,varIDX,1)=nanmean(zData);
    yStdErr(:,varIDX,1)=nanstd(zData)/sqrt(size(zData,1));
    zData=data.(VarList{varIDX}).DPrime.First3.DecZ;
    y(:,varIDX,2)=nanmean(zData);
    yStdErr(:,varIDX,2)=nanstd(zData)/sqrt(size(zData,1));
end
for condIDX=1:2
    subplot(2,1,condIDX);
    plot(data.xA,y(:,:,condIDX),'linewidth',2.5);
    hold on
    %plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','.','markersize',3,'linestyle','none');
    plot(downsample(data.xA,1),downsample(y(:,:,condIDX)-yStdErr(:,:,condIDX),1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    plot(downsample(data.xA,1),downsample(y(:,:,condIDX)+yStdErr(:,:,condIDX),1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    %plot(data.xA,y+yStdErr,'linestyle',':');
    title(['Overall zscore']);
    axis([-3 13 -2 2]);
    legend (VarList{:});
end
%% Dprime Output for Stats
VarList={'LL04','RL04','LL10','RL10'};
TimePointList={'E' 'E' 'L' 'L'};
IDArrayInc=[data.LL04.DPrime.First3.IncIDs;data.RL04.DPrime.First3.IncIDs;...
    data.LL10.DPrime.First3.IncIDs+51;data.RL10.DPrime.First3.IncIDs+51];
IDArrayDec=[data.LL04.DPrime.First3.DecIDs;data.RL04.DPrime.First3.DecIDs;...
    data.LL10.DPrime.First3.DecIDs+51;data.RL10.DPrime.First3.DecIDs+51];
TimePointArrayInc=cell(length(IDArrayInc),1);
TimePointArrayDec=cell(length(IDArrayDec),1);
DataArrayInc=NaN(length(IDArrayInc),length(data.xA));
DataArrayDec=NaN(length(IDArrayDec),length(data.xA));
xInc=1;
xDec=1;
for varIDX=1:4
    tempInc=length(data.(VarList{varIDX}).DPrime.First3.IncIDs);
    tempDec=length(data.(VarList{varIDX}).DPrime.First3.DecIDs);
    TimePointArrayInc(xInc:xInc+tempInc-1)={TimePointList{varIDX}};
    TimePointArrayDec(xDec:xDec+tempDec-1)={TimePointList{varIDX}};
    DataArrayInc(xInc:xInc+tempInc-1,:)=data.(VarList{varIDX}).DPrime.First3.IncZ;
    DataArrayDec(xDec:xDec+tempDec-1,:)=data.(VarList{varIDX}).DPrime.First3.DecZ;
    xInc=xInc+tempInc;
    xDec=xDec+tempDec;
end
%% New! Dprime Output for Stats
VarList={'LL04','RL04','LL10','RL10'};
TimePointList={'E' 'E' 'L' 'L'};
ScheduleList={'RI' 'RR' 'RI' 'RR'}
DPrime.IDArrayInc=[data.LL04.DPrime.First3.IncIDs;data.RL04.DPrime.First3.IncIDs;...
    data.LL10.DPrime.First3.IncIDs+51;data.RL10.DPrime.First3.IncIDs+51];
DPrime.IDArrayDec=[data.LL04.DPrime.First3.DecIDs;data.RL04.DPrime.First3.DecIDs;...
    data.LL10.DPrime.First3.DecIDs+51;data.RL10.DPrime.First3.DecIDs+51];
DPrime.TimePointArrayInc=cell(length(DPrime.IDArrayInc),1);
DPrime.TimePointArrayDec=cell(length(DPrime.IDArrayDec),1);
DPrime.ScheduleArrayInc=cell(length(DPrime.IDArrayInc),1);
DPrime.ScheduleArrayDec=cell(length(DPrime.IDArrayDec),1);
DPrime.DataArrayInc=NaN(length(DPrime.IDArrayInc),length(data.xA));
DPrime.DataArrayDec=NaN(length(DPrime.IDArrayDec),length(data.xA));
length1=length(data.LL04.DPrime.First3.CIDPIdx);
length2=length(data.LL10.DPrime.First3.CIDPIdx);
DPrime.Category.IDArray=[1:length1 length1+1:length1+length2]';
DPrime.Category.TimePoint=cell(length1+length2,1);
DPrime.Category.Cat1=cell(length1+length2,1);
DPrime.Category.Cat2=cell(length1+length2,1);
DPrime.Category.TimePoint(1:length1) = {'E'};
DPrime.Category.TimePoint(length1+1:length1+length2) = {'L'};
xInc=1;
xDec=1;
for varIDX=1:4
    if varIDX<=2
        cellOffset=0;
    else
        cellOffset=length(data.LL04.DPrime.First3.CIDPIdx);
    end
    if varIDX==2|| varIDX==4
        CatIDX='Cat2';
    else
        CatIDX='Cat1';
    end
    DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.IncIDs+cellOffset)={'I'};
    DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.DecIDs+cellOffset)={'D'};
    DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.NCIDs+cellOffset)={'N'};
    
    
    tempInc=length(data.(VarList{varIDX}).DPrime.First3.IncIDs);
    tempDec=length(data.(VarList{varIDX}).DPrime.First3.DecIDs);
    DPrime.tempInc=length(data.(VarList{varIDX}).DPrime.First3.IncIDs);
    DPrime.tempDec=length(data.(VarList{varIDX}).DPrime.First3.DecIDs);
    DPrime.TimePointArrayInc(xInc:xInc+tempInc-1)={TimePointList{varIDX}};
    DPrime.TimePointArrayDec(xDec:xDec+tempDec-1)={TimePointList{varIDX}};
    DPrime.ScheduleArrayInc(xInc:xInc+tempInc-1)={ScheduleList{varIDX}};
    DPrime.ScheduleArrayDec(xDec:xDec+tempDec-1)={ScheduleList{varIDX}};
    DPrime.DataArrayInc(xInc:xInc+tempInc-1,:)=data.(VarList{varIDX}).DPrime.First3.IncZ;
    DPrime.DataArrayDec(xDec:xDec+tempDec-1,:)=data.(VarList{varIDX}).DPrime.First3.DecZ;
    xInc=xInc+tempInc;
    xDec=xDec+tempDec;
end

%% DPrime Cell overlappage
overlap(1,1)=length(intersect(data.LL04.DPrime.First3.IncIDs,data.RL04.DPrime.First3.IncIDs))/length(data.LL04.DPrime.First3.IncIDs);
overlap(2,1)=length(intersect(data.LL04.DPrime.First3.IncIDs,data.RL04.DPrime.First3.IncIDs))/length(data.RL04.DPrime.First3.IncIDs);
overlap(3,1)=length(intersect(data.LL10.DPrime.First3.IncIDs,data.RL10.DPrime.First3.IncIDs))/length(data.LL10.DPrime.First3.IncIDs);
overlap(4,1)=length(intersect(data.LL10.DPrime.First3.IncIDs,data.RL10.DPrime.First3.IncIDs))/length(data.RL10.DPrime.First3.IncIDs);

overlap(1,2)=length(intersect(data.LL04.DPrime.First3.DecIDs,data.RL04.DPrime.First3.DecIDs))/length(data.LL04.DPrime.First3.DecIDs);
overlap(2,2)=length(intersect(data.LL04.DPrime.First3.DecIDs,data.RL04.DPrime.First3.DecIDs))/length(data.RL04.DPrime.First3.DecIDs);
overlap(3,2)=length(intersect(data.LL10.DPrime.First3.DecIDs,data.RL10.DPrime.First3.DecIDs))/length(data.LL10.DPrime.First3.DecIDs);
overlap(4,2)=length(intersect(data.LL10.DPrime.First3.DecIDs,data.RL10.DPrime.First3.DecIDs))/length(data.RL10.DPrime.First3.DecIDs);
%% DPrime Cell overlappage reversed
overlap(1,1)=length(intersect(data.LL04.DPrime.First3.IncIDs,data.RL04.DPrime.First3.DecIDs))/length(data.LL04.DPrime.First3.IncIDs);
overlap(2,1)=length(intersect(data.LL04.DPrime.First3.DecIDs,data.RL04.DPrime.First3.IncIDs))/length(data.RL04.DPrime.First3.IncIDs);
overlap(3,1)=length(intersect(data.LL10.DPrime.First3.IncIDs,data.RL10.DPrime.First3.DecIDs))/length(data.LL10.DPrime.First3.IncIDs);
overlap(4,1)=length(intersect(data.LL10.DPrime.First3.DecIDs,data.RL10.DPrime.First3.IncIDs))/length(data.RL10.DPrime.First3.IncIDs);

overlap(1,2)=length(intersect(data.LL04.DPrime.First3.DecIDs,data.RL04.DPrime.First3.IncIDs))/length(data.LL04.DPrime.First3.DecIDs);
overlap(2,2)=length(intersect(data.LL04.DPrime.First3.IncIDs,data.RL04.DPrime.First3.DecIDs))/length(data.RL04.DPrime.First3.DecIDs);
overlap(3,2)=length(intersect(data.LL10.DPrime.First3.DecIDs,data.RL10.DPrime.First3.IncIDs))/length(data.LL10.DPrime.First3.DecIDs);
overlap(4,2)=length(intersect(data.LL10.DPrime.First3.IncIDs,data.RL10.DPrime.First3.DecIDs))/length(data.RL10.DPrime.First3.DecIDs);
%% DPrime Cell overlappage - either direction
LL04All=[data.LL04.DPrime.First3.IncIDs;data.LL04.DPrime.First3.DecIDs];
RL04All=[data.RL04.DPrime.First3.IncIDs;data.RL04.DPrime.First3.DecIDs];
LL10All=[data.LL10.DPrime.First3.IncIDs;data.LL10.DPrime.First3.DecIDs];
RL10All=[data.RL10.DPrime.First3.IncIDs;data.RL10.DPrime.First3.DecIDs];
overlapAll(1,1)=length(intersect(LL04All,RL04All))/length(LL04All);
overlapAll(2,1)=length(intersect(LL04All,RL04All))/length(RL04All);
overlapAll(1,2)=length(intersect(LL10All,RL10All))/length(LL10All);
overlapAll(2,2)=length(intersect(LL10All,RL10All))/length(RL10All);
%% 
overlap(1,:)=[length(intersect(data.LL04.DPrime.First3.IncIDs,data.RL04.DPrime.First3.IncIDs)) length(data.LL04.DPrime.First3.IncIDs)];
overlap(2,:)=[length(intersect(data.LL04.DPrime.First3.IncIDs,data.RL04.DPrime.First3.IncIDs)) length(data.RL04.DPrime.First3.IncIDs)];
overlap(3,:)=[length(intersect(data.LL10.DPrime.First3.IncIDs,data.RL10.DPrime.First3.IncIDs)) length(data.LL10.DPrime.First3.IncIDs)];
overlap(4,:)=[length(intersect(data.LL10.DPrime.First3.IncIDs,data.RL10.DPrime.First3.IncIDs)) length(data.RL10.DPrime.First3.IncIDs)];

overlap(5,:)=[length(intersect(data.LL04.DPrime.First3.DecIDs,data.RL04.DPrime.First3.DecIDs)) length(data.LL04.DPrime.First3.DecIDs)];
overlap(6,:)=[length(intersect(data.LL04.DPrime.First3.DecIDs,data.RL04.DPrime.First3.DecIDs)) length(data.RL04.DPrime.First3.DecIDs)];
overlap(7,:)=[length(intersect(data.LL10.DPrime.First3.DecIDs,data.RL10.DPrime.First3.DecIDs)) length(data.LL10.DPrime.First3.DecIDs)];
overlap(8,:)=[length(intersect(data.LL10.DPrime.First3.DecIDs,data.RL10.DPrime.First3.DecIDs)) length(data.RL10.DPrime.First3.DecIDs)];

%% DPrime Cell overlappage reversed
overlapR(1,:)=[length(intersect(data.LL04.DPrime.First3.IncIDs,data.RL04.DPrime.First3.DecIDs)) length(data.LL04.DPrime.First3.IncIDs)];
overlapR(2,:)=[length(intersect(data.LL04.DPrime.First3.DecIDs,data.RL04.DPrime.First3.IncIDs)) length(data.RL04.DPrime.First3.IncIDs)];
overlapR(3,:)=[length(intersect(data.LL10.DPrime.First3.IncIDs,data.RL10.DPrime.First3.DecIDs)) length(data.LL10.DPrime.First3.IncIDs)];
overlapR(4,:)=[length(intersect(data.LL10.DPrime.First3.DecIDs,data.RL10.DPrime.First3.IncIDs)) length(data.RL10.DPrime.First3.IncIDs)];

overlapR(5,:)=[length(intersect(data.LL04.DPrime.First3.DecIDs,data.RL04.DPrime.First3.IncIDs)) length(data.LL04.DPrime.First3.DecIDs)];
overlapR(6,:)=[length(intersect(data.LL04.DPrime.First3.IncIDs,data.RL04.DPrime.First3.DecIDs)) length(data.RL04.DPrime.First3.DecIDs)];
overlapR(7,:)=[length(intersect(data.LL10.DPrime.First3.DecIDs,data.RL10.DPrime.First3.IncIDs)) length(data.LL10.DPrime.First3.DecIDs)];
overlapR(8,:)=[length(intersect(data.LL10.DPrime.First3.IncIDs,data.RL10.DPrime.First3.DecIDs)) length(data.RL10.DPrime.First3.DecIDs)];
%% DPrime Cell overlappage - either direction
LL04All=[data.LL04.DPrime.First3.IncIDs;data.LL04.DPrime.First3.DecIDs];
RL04All=[data.RL04.DPrime.First3.IncIDs;data.RL04.DPrime.First3.DecIDs];
LL10All=[data.LL10.DPrime.First3.IncIDs;data.LL10.DPrime.First3.DecIDs];
RL10All=[data.RL10.DPrime.First3.IncIDs;data.RL10.DPrime.First3.DecIDs];
overlapAll(1,:)=[length(intersect(LL04All,RL04All)) length(LL04All)];
overlapAll(2,:)=[length(intersect(LL04All,RL04All)) length(RL04All)];
overlapAll(3,:)=[length(intersect(LL10All,RL10All)) length(LL10All)];
overlapAll(4,:)=[length(intersect(LL10All,RL10All)) length(RL10All)];
%% Getting the pie charts
VarList={'LL04','RL04','LL10','RL10'};
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    pie(data.(VarList{varIDX}).DPrime.First3.pieData,{['Inc: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(1))] ...
        ['NoChange: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(2))] ...
        ['Dec: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(3))]});
end
%% Inc and Dec
VarList={'LL04','RL04','LL10','RL10'};
numCells=10;
preRange=find(data.xA>=-1.7 & data.xA<0);
postRange=find(data.xA>=0 & data.xA<=1.7);
for varIDX=1:4
    data.(VarList{varIDX}).mySort.IncZ=zscore(data.(VarList{varIDX}).MlMtx(data.(VarList{varIDX}).mySort.idx(1:numCells),:)')';
    data.(VarList{varIDX}).mySort.IncZPre1_7=data.(VarList{varIDX}).mySort.IncZ(:,preRange);
    data.(VarList{varIDX}).mySort.IncZPost1_7=data.(VarList{varIDX}).mySort.IncZ(:,postRange);
    data.(VarList{varIDX}).mySort.DecZ=zscore(data.(VarList{varIDX}).MlMtx(data.(VarList{varIDX}).mySort.idx(end-(numCells-1):end),:)')';
    data.(VarList{varIDX}).mySort.DecZPre1_7=data.(VarList{varIDX}).mySort.DecZ(:,preRange);
    data.(VarList{varIDX}).mySort.DecZPost1_7=data.(VarList{varIDX}).mySort.DecZ(:,postRange);
end
%% Inc and Dec Plot and data out
VarList={'LL04','RL04','LL10','RL10'};
for varIDX=1:4
    output.Inc(varIDX,:)=mean(data.(VarList{varIDX}).mySort.IncZ);
    output.Inc(varIDX+4,:)=std(data.(VarList{varIDX}).mySort.IncZ)/sqrt(numCells);
    output.IncPre1_7(varIDX,:)=mean(data.(VarList{varIDX}).mySort.IncZPre1_7);
    output.IncPre1_7(varIDX+4,:)=std(data.(VarList{varIDX}).mySort.IncZPre1_7)/sqrt(numCells);
    output.IncPost1_7(varIDX,:)=mean(data.(VarList{varIDX}).mySort.IncZPost1_7);
    output.IncPost1_7(varIDX+4,:)=std(data.(VarList{varIDX}).mySort.IncZPost1_7)/sqrt(numCells);
    
    output.Dec(varIDX,:)=mean(data.(VarList{varIDX}).mySort.DecZ);
    output.Dec(varIDX+4,:)=std(data.(VarList{varIDX}).mySort.DecZ)/sqrt(numCells);
    output.DecPre1_7(varIDX,:)=mean(data.(VarList{varIDX}).mySort.DecZPre1_7);
    output.DecPre1_7(varIDX+4,:)=std(data.(VarList{varIDX}).mySort.DecZPre1_7)/sqrt(numCells);
    output.DecPost1_7(varIDX,:)=mean(data.(VarList{varIDX}).mySort.DecZPost1_7);
    output.DecPost1_7(varIDX+4,:)=std(data.(VarList{varIDX}).mySort.DecZPost1_7)/sqrt(numCells);
end
%% lines
figure;
CondList={'Inc','Dec'};
meanRange=1:4;
stderrRange=5:8;
for condIDX=1:2
    subplot(2,1,condIDX);
    y=(output.(CondList{condIDX})(meanRange,:));
    yStdErr=(output.(CondList{condIDX})(stderrRange,:));
    plot(data.xA,y,'linewidth',2.5);
    hold on
    %plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','.','markersize',3,'linestyle','none');
    plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    plot(downsample(data.xA,1),downsample(y+yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    %plot(data.xA,y+yStdErr,'linestyle',':');
    title(['Overall zscore: ' CondList{condIDX}]);
    axis([-3 3 -1.5 2.5]);
    legend (VarList{:});
end

%% modify to get bars
CondList={'IncPre1_7','IncPost1_7','DecPre1_7','DecPost1_7'};
meanRange=1:4;
stderrRange=5:8;
for varIDX=1:4
    subplot(2,2,varIDX);
    y=(output.(CondList{varIDX})(meanRange,:))';
    yStdErr=(output.(CondList{varIDX})(stderrRange,:))';
    plot(data.xA,y,'linewidth',.8);
    hold on
    %plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','.','markersize',3,'linestyle','none');
    plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    plot(downsample(data.xA,1),downsample(y+yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
    %plot(data.xA,y+yStdErr,'linestyle',':');
    title(['Overall zscore: ' CondList{varIDX}]);
    axis([-3 10 -1 1])
end
