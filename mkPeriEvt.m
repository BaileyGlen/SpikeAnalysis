function [eventStruct]=mkPeriEvt ...
    (evtTrigger,DirList,preEvt, postEvt, rasterBin)
% MKPERIEVT generate the peri-event histogram around an event.
% evtTrigger - The event to look for in the dataset, plus the type of
% analysis.
% clr        - vestigial
% DirList    - the list of files to process with this evt trigger
% zScore     - vestigial
% PLOT       - vestigial but is definitely worth reimplementing
% preEvt     - time before event to capture in Secs
% postEvt    - time after event to capture
% rasterBin  - the size of the individual bin
%Event Triggers:
%

%left lever  = RI  or  LL
%right elver = RR  or  RL

% clear all
% close all
% PLOT=1;
% % evtTrigger={'LL_R';'RL_R'}; clr='k';
% % evtTrigger={'RL_R'}; clr='r';
% evtTrigger={'LL_R'}; clr='b';
% selDay=04;
% zScore='zNo';
% rasterBin=0.05;
% preEvt=3;            % time prior to Event in sec
% postEvt=10;           % time post Event in sec

%         %DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m22D04_03';'m25D04_03.mat'};
%         DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m21D04_03';'m22D04_03'};
%     case 10
%         DirList={'m13D04_03.mat';'m14D10_03.mat';'m25D10_03.mat'};
%     case 410
%         DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m21D04_03.mat';'m22D04_03.mat';'m25D04_03.mat';...
%             'm13D10_03.mat';'m14D04_03.mat';'m21D10_03.mat';'m22D10_03.mat';'m25D10_03.mat'};


%%
consumptiondelay = 15;

eventStruct.AnimalIDList={};
eventStruct.AnimalIDList_Dropped = {};
numFiles = length(find(~cellfun(@isempty,DirList)));
for XX=1:numFiles
    %% Load data
    load(DirList{XX});
    data.Freq = 40000;
    if ~exist('data','var')
        data=struct();
        data.SessionLength = 30;
        data.numLevers=2;
    elseif isfield (data,'Lapish')
        behaveEvt_Raw = data.Lapish.behaveEvt_Raw;
        behaveEvtTm_Raw = data.Lapish.behaveEvtTm_Raw;
        STMtx = data.Lapish.STMtx;
        mask = data.goodDetectionMask;
        %set the first and last .2s to bad to remove filtfilt ripples
        mask(1:400*.2) = 0;
        mask(end-(400*.1):end) = 0;
        %find the transitions to 0, and remove the points before. then do
        %to 1
        zerovals = find(diff(mask) == -1)+1;
        for x = 1:length(zerovals)
            mask(zerovals(x)-(.05*400):zerovals(x)) = 0;
        end
        %% these onevals aren't being used below and that makes no sense
        onevals = find(diff(mask) == +1)+1;
        for x = 1:length(onevals)
            mask(onevals(x):onevals(x)+(.05*400)) = 0;
        end
        maskup = upsample(mask,100);
        for x = 1:data.SessionLength*60*data.Freq/100
            maskup(x*100-99:x*100) = mask (x);
        end
        mask = maskup;
    end
    data.AnimalID = ['M' DirList{XX}(2:3)];
    
    SurSpike={};
    % Mean Firing Rates
    %% Consumption
    if strcmp(evtTrigger,'RL_C')
        Press_raw=strmatch('RL_R',behaveEvt_Raw);
        if isempty(Press_raw)
            Press_raw=strcmp('RR_R',behaveEvt_Raw);
            display('Had to correct one RL to RR');
        end
        Press_ts=behaveEvtTm_Raw(Press_raw);
        BB_raw=strcmp('BB',behaveEvt_Raw);
        BB_ts=behaveEvtTm_Raw(BB_raw);
        k=unique(cell2mat(arrayfun (@(x) find(BB_ts>x & BB_ts<x + consumptiondelay,1), Press_ts, 'UniformOutput', false)));
        hld=BB_ts(k);
    elseif strcmp(evtTrigger,'LL_C')
        Press_raw=strcmp('LL_R',behaveEvt_Raw);
        Press_ts=behaveEvtTm_Raw(Press_raw);
        BB_raw=strcmp('BB',behaveEvt_Raw);
        BB_ts=behaveEvtTm_Raw(BB_raw);
        k=unique(cell2mat(arrayfun (@(x) find(BB_ts>x & BB_ts<x + consumptiondelay,1), Press_ts, 'UniformOutput', false)));
        hld=BB_ts(k);
    elseif strcmp(evtTrigger(4:end),'Sucrose_C')
        Press_raw=strcmp('RF',behaveEvt_Raw);
        %         if isempty(Press_raw)
        %             Press_raw=strcmp('RR_R',behaveEvt_Raw);
        %             display('Had to correct one RL to RR');
        %         end
        Press_ts=behaveEvtTm_Raw(Press_raw);
        BB_raw=strcmp('BB',behaveEvt_Raw);
        BB_ts=behaveEvtTm_Raw(BB_raw);
        k=unique(cell2mat(arrayfun (@(x) find(BB_ts>x & BB_ts<x+ consumptiondelay,1), Press_ts, 'UniformOutput', false)));
        hld=BB_ts(k);
        %% Isolated Presses
    elseif strcmp(evtTrigger,'RL_I')
        allEvents = {'BB' 'RF' 'RL_R' 'LL_R' 'RL_U' 'LL_U'};
        Press_raw=strmatch('RL_U',behaveEvt_Raw);
        if isempty(Press_raw)
            display('Had to correct one RL to RR');
            Press_raw=strcmp('RR_U',behaveEvt_Raw);
            allEvents = {'BB' 'RF' 'RR_R' 'LL_R' 'RR_U' 'LL_U'};
        end
        Press_ts=behaveEvtTm_Raw(Press_raw);
        Other_raw= cellfun(@(y) max(y), cellfun(@(x) strcmp(allEvents,x),behaveEvt_Raw,'UniformOutput', false) );
        Other_ts=behaveEvtTm_Raw(Other_raw);
        k=cellfun(@(y) isempty(y), arrayfun (@(x) find(Other_ts>(x-3) & Other_ts<(x+3) & abs(Other_ts-x)>.0001), Press_ts, 'UniformOutput', false)  );
        hld=Press_ts(k);
        postEvt=3;
    elseif strcmp(evtTrigger,'LL_I')
        allEvents = {'BB' 'RF' 'RL_R' 'LL_R' 'RL_U' 'LL_U'};
        Press_raw=strcmp('LL_U',behaveEvt_Raw);
        Press_ts=behaveEvtTm_Raw(Press_raw);
        Other_raw= cellfun(@(y) max(y), cellfun(@(x) strcmp(allEvents,x),behaveEvt_Raw,'UniformOutput', false) );
        Other_ts=behaveEvtTm_Raw(Other_raw);
        k=cellfun(@(y) isempty(y), arrayfun (@(x) find(Other_ts>(x-3) & Other_ts<(x+3) & abs(Other_ts-x)>.0001), Press_ts, 'UniformOutput', false)  );
        hld=Press_ts(k);
        postEvt=3;
        %     elseif strcmp(evtTrigger,'RL_Sucrose_I')
        %         allEvents = {'BB' 'RF' 'RL_U' 'LL_U'};
        %         Press_raw=strmatch('RL_U',behaveEvt_Raw);
        %         if isempty(Press_raw)
        %             display('Had to correct one RL to RR');
        %             Press_raw=strcmp('RR_U',behaveEvt_Raw);
        %             allEvents = {'BB' 'RF' 'RR_U' 'LL_U'};
        %         end
        %         Press_ts=behaveEvtTm_Raw(Press_raw);
        %         Other_raw= cellfun(@(y) max(y), cellfun(@(x) strcmp(allEvents,x),behaveEvt_Raw,'UniformOutput', false) );
        %         Other_ts=behaveEvtTm_Raw(Other_raw);
        %         k=cellfun(@(y) isempty(y), arrayfun (@(x) find(Other_ts>(x-3) & Other_ts<(x+3) & abs(Other_ts-x)>.0001), Press_ts, 'UniformOutput', false)  );
        %         hld=Press_ts(k);
        %         postEvt=3;
        %     elseif strcmp(evtTrigger,'LL_Sucrose_I')
        %         allEvents = {'BB' 'RF' 'RL_U' 'LL_U'};
        %         Press_raw=strmatch('LL_U',behaveEvt_Raw);
        %         Press_ts=behaveEvtTm_Raw(Press_raw);
        %         Other_raw= cellfun(@(y) max(y), cellfun(@(x) strcmp(allEvents,x),behaveEvt_Raw,'UniformOutput', false) );
        %         Other_ts=behaveEvtTm_Raw(Other_raw);
        %         k=cellfun(@(y) isempty(y), arrayfun (@(x) find(Other_ts>(x-3) & Other_ts<(x+3) & abs(Other_ts-x)>.0001), Press_ts, 'UniformOutput', false)  );
        %         hld=Press_ts(k);
        %         postEvt=3;
        %% Reinforcer Delivery - rr= random ratio
    elseif strcmp(evtTrigger(4:end),'Sucrose')
        k=strcmp('RF',behaveEvt_Raw);
        hld=behaveEvtTm_Raw(k);
    elseif strcmp(evtTrigger(4), 'R')     % assumes 'R' for reinforcer
        k=strmatch(evtTrigger,behaveEvt_Raw);
        if isempty(k)
            if strcmp(evtTrigger(1:2),'RL')
                origTrigger = evtTrigger;
                evtTrigger(1:2) = 'RR';
                display('Had to correct one RL to RR');
                k=strcmp(evtTrigger,behaveEvt_Raw);
            end
        end
        hld=behaveEvtTm_Raw(k);
        if strcmp(evtTrigger(4), 'R')
            k2=strcmp([evtTrigger(1:3) 'U'], behaveEvt_Raw);
            hld2=behaveEvtTm_Raw(k2);
        end
        if exist('origTrigger') 
            evtTrigger = origTrigger;
        end
    else
        error(['The evtTrigger provided is inaccurate:,' evtTrigger]);
    end
    %selEvt=sort(cell2mat(hld'));
    selEvt=hld;
    if (exist('k2'))
        selEvt2=hld2;
    end
    % Tacking on behavrioal events to STMtx in last column to build rasters
    STMtx(:,end+1)=repmat(nan,size(STMtx,1),1);
    STMtx(1:length(selEvt),end)=selEvt;
    
    
    
    % Building Raster (with selected beahavioral events in last column)
    rasterTm=[0:rasterBin:max(STMtx(:))];
    hRaster=histc(STMtx,rasterTm);
    k=find(hRaster(:,end)~=0);
    % Making histogram of artifacts to reject
    artHist=histc(find(mask==0)/40000,rasterTm);
    if ~isempty(k) %k isn't empty
        
        for i=1:size(hRaster,2)-1;   % Number of spike trains
            for j=1:length(k);     % Number of events
                if k(j)-preEvt*(1/rasterBin)>0 &&  k(j)+postEvt*(1/rasterBin)<length(hRaster);
                    pEvt{XX,1}{i,1}(j,:)=smooth(hRaster(k(j)-preEvt*(1/rasterBin):k(j)+postEvt*(1/rasterBin),i),5);
                    if i==1 && (exist('k2'))
                        SurSpike{XX,1}(j,1)=length(find(selEvt2>=selEvt(j)-3 & selEvt2<selEvt(j)));
                        SurSpike{XX,1}(j,2)=length(find(selEvt2>selEvt(j) & selEvt2<=selEvt(j)+3));
                    end
                    %%% BG - If histc bin has 1, set pEvt to NaN - this could be just
                    %%% the bin, or the entire row
                    needNaN = find(artHist(k(j)-preEvt*(1/rasterBin):k(j)+postEvt*(1/rasterBin)));
                    %set the misisng Bins to NaN
                    %arrayfun(@(QQ) pEvt{QQ,1}{i,1}(j,i) = NaN, find(artHist(k(j)-preEvt*(1/rasterBin):k(j)+postEvt*(1/rasterBin))));
                    %                 for x= needNaN
                    %                     pEvt{XX,1}{i,1}(j,needNaN)= NaN;
                    %                 end
                    % set the entire Event to NaN
                    if needNaN
                        pEvt{XX,1}{i,1}(j,:)= NaN;
                    end
                    %%% BG - end
                end;
            end;
            %buildBaseline pEvt
            %[baselineRange]=BaselineGen(behaveEvtTm_Raw, mask);
            %behaveEvtTm_Raw(end+1,1) = 900;
            %behaveEvtTm_Raw = sort(behaveEvtTm_Raw);
            if mean(k<data.SessionLength*60)
                diffArray = diff(behaveEvtTm_Raw(behaveEvtTm_Raw<(data.SessionLength*60/2)));
            else
                diffArray = diff(behaveEvtTm_Raw(behaveEvtTm_Raw>(data.SessionLength*60/2)));
            end
            emptyEpochsIDX = find(diffArray>20);
            for j=1:size(emptyEpochsIDX,1)
                centerpt = round(behaveEvtTm_Raw(emptyEpochsIDX(j))+ diffArray(emptyEpochsIDX(j))/2);
                pEvt_base{XX,1}{i,1}(j,:)=smooth(hRaster((centerpt-6.5)*(1/rasterBin):(centerpt+6.5)*(1/rasterBin),i),5);
                needNaN = find(artHist((centerpt-6.5)*(1/rasterBin):(centerpt+6.5)*(1/rasterBin)));
                %set the misisng Bins to NaN
                %arrayfun(@(QQ) pEvt{QQ,1}{i,1}(j,i) = NaN, find(artHist(k(j)-preEvt*(1/rasterBin):k(j)+postEvt*(1/rasterBin))));
                for x= needNaN
                    pEvt_base{XX,1}{i,1}(j,needNaN)= NaN;
                end
            end
        end
        if all(isnan(pEvt{XX,1}{1}))
            pEvt{XX,1} = [];
            pEvt_base{XX,1} = [];
            eventStruct.AnimalIDList_Dropped{length(eventStruct.AnimalIDList_Dropped)+1}={data.AnimalID};
        else
            eventStruct.AnimalIDList{1,XX} = data.AnimalID;
        end
        
    else
        pEvt{XX,1} = [];
        pEvt_base{XX,1} = [];
        eventStruct.AnimalIDList_Dropped{length(eventStruct.AnimalIDList_Dropped)+1}={data.AnimalID};
    end
end
if ~isempty(pEvt)
    for i=1:numFiles;
        dMtx{i,1}=cell2mat(pEvt{i});
        for j=1:length(pEvt{i})
            MndMtx{i,1}(j,:)=nanmean(pEvt{i}{j},1);
        end
        for j=1:length(pEvt_base{i})
            BldMtx{i,1}(j,:)=nanmean(pEvt_base{i}{j},1);
        end
    end
    eventStruct.Mtx=cell2mat(dMtx);
    eventStruct.MlMtx = cell2mat(MndMtx);
    eventStruct.BlMtx=cell2mat(BldMtx);
    
end
eventStruct.pEvt = pEvt;
eventStruct.pEvt_base = pEvt_base;
eventStruct.SurSpike = SurSpike;

xA=[-1*preEvt:rasterBin:postEvt];

% %% Plotting tool
% if PLOT==1;
%     figure(1)
%     switch zScore
%         case 'zYes'
%             plot(xA,nanmean(zscore(Mtx')'),clr);hold on;
%         case 'zNo'
%             plot(xA,nanmean(Mtx),clr);hold on;
%     end
%
%     figure();
%     switch zScore
%         case 'zYes'
%             imagesc(xA,1:length(Mtx),zscore(Mtx')',[-2 2]);
%         case 'zNo'
%             imagesc(xA,1:size(Mtx,2),Mtx,[0 3]);
%     end;
%     colormap('jet');
%     figure(3);
%     switch zScore
%         case 'zYes'
%             plot(xA,nanmean(zscore(MnMtx')'),clr);hold on;
%         case 'zNo'
%             plot(xA,nanmean(MnMtx),clr);hold on;
%     end
%
%     figure();
%     switch zScore
%         case 'zYes'
%             imagesc(xA,1:length(MnMtx),zscore(MnMtx')',[-2 2]);
%         case 'zNo'
%             imagesc(xA,1:size(MnMtx,2),MnMtx,[0 3]);
%     end;
%     colormap('jet');
% end;
