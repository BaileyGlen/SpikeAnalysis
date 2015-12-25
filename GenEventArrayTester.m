var1={'LL_R' 'RL_R'};
var2={'LL_R' 'RL_R' 'BB'  'LL_U' 'RL_U'};
eventCellArray=cell(2,1);
eventCellArray{1}=[1; 2; 3; 10; 20; ...
    30; 31; 32; 33; 34; ...
    35; 50; 55; 60; 65];
eventCellArray{2}={'RL_U'; 'BB'; 'LL_U'; 'LL_R'; 'BB'; ...
    'LL_U';  'LL_R';  'BB';  'RL_R';  'LL_R'; ...
    'BB';  'BB';  'BB';  'LL_R'; 'LL_R'};
GenEventArray(eventCellArray,var1,var2,'opt','iso');
%%
var1={'LL_R'};
var2={'BB'};
eventCellArray=cell(2,1);
eventCellArray{1}=[1; 2; 3; 10; 20; ...
    30; 31; 32; 33; 34; ...
    35; 50; 55; 60; 65];
eventCellArray{2}={'LL_R'; 'BB'; 'LL_U'; 'LL_R'; 'BB'; ...
    'LL_U';  'LL_R';  'LL_R';  'RL_R';  'LL_U'; ...
    'BB';  'BB';  'BB';  'LL_R'; 'LL_R'};
[~,stats]=GenEventArray(eventCellArray,var1,var2,'opt','firstAfter');

%%
% mnSpk = 0;
% mxSpk=1000;
% rasterBin=.1;      % Size of bines for raster
% preEvt=3;            % time prior to Event in sec
% postEvt=3;           % time post Event in sec
% PLOT=0; 
% nmPC=4;
%tempRange = [1:((preEvt+postEvt)/rasterBin)+1];
DayVar={4,4,10,10};
%EventVar={'LL_I','RL_I','LL_I','RL_I'};
EventVar={'LL_R','RL_R','LL_R','RL_R'};
%EventVar={'LL_C','RL_C','LL_C','RL_C'};
StructVar={'LL04','RL04','LL10','RL10'};
%% Create the data Structure - Reinforced
%dataStruct.xA=[-1*preEvt:rasterBin:postEvt];
for x=1:4
    Day=DayVar{x};
    %bEvt = EventVar{x};
    switch Day
        case 4
            %DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m22D04_03';'m25D04_03.mat'};
            %DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m21D04_03';'m22D04_03';'m25D04_03.mat'};%All
            DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m22D04_03';'m25D04_03.mat'}; %RfCon
        case 10
            %DirList={'m13D10_03.mat';'m14D10_03.mat';'m16D10_03.mat';'m21D10_03';'m22D10_03';'m25D10_03.mat'};%All
            %DirList={'m13D10_03.mat';'m14D10_03.mat';'m21D10_03';'m22D10_03';'m25D10_03.mat'};%RfIso
            DirList={'m13D10_03.mat';'m14D10_03.mat';'m21D10_03.mat';'m22D10_03';'m25D10_03.mat'};%RfCon
        case 410
            DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m21D04_03.mat';'m22D04_03.mat';'m25D04_03.mat';...
                'm13D10_03.mat';'m14D04_03.mat';'m21D10_03.mat';'m22D10_03.mat';'m25D10_03.mat'};
    end
    
    for XX=1:length(DirList)
        % Load data
        load(DirList{XX});
        [~,stats]=GenEventArray({behaveEvtTm_Raw; behaveEvt_Raw},EventVar{x},'BB','opt','firstAfter');
        data.(StructVar{x}).Rf2BB.Lag{XX,1}=stats.Lag;
        %dataStruct.(StructVar{x})=getPCA(dataStruct.(StructVar{x}));
    end
end
%%
VarList={'LL04','RL04','LL10','RL10'};
for varIDX=1:4
    totCount=1;
    for CellIDX=1:length(data.(VarList{varIDX}).Rf2BB.Lag)
        multVar=length(data.(VarList{varIDX}).pEvt{CellIDX});
        temp(totCount:totCount+multVar,varIDX)=nanmean(data.(VarList{varIDX}).Rf2BB.Lag{CellIDX}(~isnan(data.(VarList{varIDX}).pEvt{CellIDX}{1}(:,1))));
        totCount = totCount + multVar;
    end
end
%%
VarList={'LL04','RL04','LL10','RL10'};
temp=[];
totCount=1;
for varIDX=1:4
    for CellIDX=1:length(data.(VarList{varIDX}).Rf2BB.Lag)
        multVar=length(data.(VarList{varIDX}).pEvt{CellIDX});%%%no pevt!!!!!!
        numTrials=length(data.(VarList{varIDX}).Rf2BB.Lag{CellIDX});
        curRange=totCount:totCount+numTrials-1;
        temp(curRange,1)=varIDX;%CondID
        temp(curRange,2)=CellIDX;%AnimalID
        temp(curRange,3)=multVar;%CellCount
        temp(curRange,4)=curRange-totCount+1;%TrialID
        %Rf2BB
        temp(curRange,5)= [data.(VarList{varIDX}).Rf2BB.Lag{CellIDX}];
        %0=Bad 1=Good
        temp(curRange,6)=~isnan(data.(VarList{varIDX}).pEvt{CellIDX}{1}(:,1));
        
        totCount = totCount + numTrials;
        
    end
end