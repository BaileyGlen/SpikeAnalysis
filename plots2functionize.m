VarList={'LL04','RL04','LL10','RL10'};

varIDX =4;
PC = 1;
figure;
[sorted,idx] = sort(data.(VarList{varIDX}).PCA.score(:,PC));
imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);

%% this is the basic 4 way PCA plot. this could be adapted to be useful for lots of things.
VarList={'LL04','RL04','LL10','RL10'};
PC = 4;
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    [sorted,idx] = sort(data.(VarList{varIDX}).PCA.score(:,PC));
    imagesc(data.xA,1:size(sorted,1),(zscore(data.(VarList{varIDX}).MlMtx(idx,:)')'),[-2 2]);
    title(['PC: ' num2str(PC) '   ' VarList{varIDX}]);
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
    axis([-3 10 -1 1])
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
    axis([-3 10 -1 1]);
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
%% Getting the pie chart data from the dprime
VarList={'LL04','RL04','LL10','RL10'};
for varIDX=1:4
    tempZ=zscore(data.(VarList{varIDX}).MlMtx')';
    SIG=data.(VarList{varIDX}).DPrime.First3.SIG;
    CIDPIdx=data.(VarList{varIDX}).DPrime.First3.CIDPIdx;
    newIDs=sort(CIDPIdx(SIG));
    %     tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    %     [sorted, idx] = sort(mean( tempZ(:,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(:,data.xA>=0 & data.xA<=3),2 ));
    
    tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    [sorted, idx] = sort(mean( tempZ(newIDs,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(newIDs,data.xA>=0 & data.xA<=3),2 ));
    data.(VarList{varIDX}).DPrime.First3.pieData(1,1)=length(find(sorted<0));
    data.(VarList{varIDX}).DPrime.First3.pieData(1,3)=length(find(sorted>0));
    data.(VarList{varIDX}).DPrime.First3.pieData(1,2)=size(tempZ,1)-size(SIG,1);
end

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
