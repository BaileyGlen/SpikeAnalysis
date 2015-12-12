

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
%% this is the overall means, zscored
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    data2plot = nanmean(zscore(data.(VarList{varIDX}).MlMtx')');
    plot(data.xA,data2plot);
    title(['Overall zscore: ' VarList{varIDX}]);
    axis([-3 10 -1 1])
end

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