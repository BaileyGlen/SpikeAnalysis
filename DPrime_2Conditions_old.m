function data = DPrime_2Conditions(data,type)%% NEWEST VERSION OF SELECTIVITY_ D.N.Linsenbardt (November 2015) dlinsen1@gmail.com
close all;
for t=1:2; % if this were like, data.CondList{}, it could read in the correct number automatically
    clearvars ('-except', 'data', 't','type');
    %cd 'D:\Users\Bailey\Documents\Dropbox\Mouse MEA\Mouse MEA\Mat\CL';
    if t==1;
        %load('LLD04_pEvt02.mat');
        curVar = 'LL';
    elseif t==2;
        %load('LLD10_pEvt02.mat');
        curVar = 'RL';
    end;
    
    % GETTING INTO SAME FORMAT AS MY PEVTMTX
    pEvt=data.(curVar).pEvt';
    
    for i=1:size(pEvt,2);
        NEWpEvt{i}=pEvt{i}';
    end;
    
    pEvtMtx=NEWpEvt;
    
    %% Assign DPrime Ranges based on type
    if type=='C' %consumption
        %these are reversed pre vs post in order to allow inc and dec to stay
        %the same
        data.preRange=data.xA>=8 & data.xA<=13;
        data.postRange=data.xA>=0 & data.xA<=5;
    elseif (type=='P' || type=='I') %Press
        data.preRange=data.xA>=-3 & data.xA<=0;
        data.postRange=data.xA>=0 & data.xA<=3;
    else error('myApp:argChk','Invalid Input for "Type"')
    end
    
    %% Initialize struct for dataset
    numRows = size(data.(curVar).MlMtx,1);
    struct4dataset.SessionType = cell(numRows,1);
    struct4dataset.Schedule = cell(numRows,1);
    struct4dataset.AnimalID = cell(numRows,1);
    struct4dataset.CellID = cell(numRows,1);
    struct4dataset.numEvents = nan(numRows,1);
    struct4dataset.numAnalyzedEvents = nan(numRows,1);
    struct4dataset.zScorePre = nan(numRows,1);
    struct4dataset.zScorePost = nan(numRows,1);
    struct4dataset.DPrimeSig = nan(numRows,1);
    struct4dataset.DPrimeCategory = cell(numRows,1);
    struct4dataset.zScoreOverall = nan(numRows,length(data.xA));
    
    %% Fill In some dataset values
    struct4dataset.SessionType(:) = {data.SessionType};
    struct4dataset.Schedule(:) = {curVar};
    %% Get the DPrimes
    m=1;
    h = waitbar(0,'Initializing waitbar...');  %
    rowCount = 1;
    for i=1:size(pEvtMtx,2); %      
        PercentComplete = i/6;  waitbar(PercentComplete,h,sprintf('%f%% Complete...D PRIME + ROC',PercentComplete));  %

        for j=1:size(pEvtMtx{i},2) %
            
            % Assign Animal and Cell ID
            CellId{m} = [int2str(i) '_' int2str(j)];m=m+1;
            struct4dataset.AnimalID{rowCount} = data.(curVar).AnimalIDList{i};
            struct4dataset.CellID{rowCount} = j;
            struct4dataset.numEvents(rowCount,1) = size(pEvtMtx{i}{j},1);
            
            struct4dataset.numAnalyzedEvents(rowCount,1) = length(find((~isnan(pEvt{i}{1}(:,1)))));  
            
            
            % Overall mean, and Individual Means
            Cellmean=zscore(nanmean(pEvtMtx{i}{j},1)); %
            CellmeanPre=nanmean(Cellmean(:,data.preRange),2); %
            CellmeanPost=nanmean(Cellmean(:,data.postRange),2); %
            struct4dataset.zScoreOverall(rowCount,:) = Cellmean;
            struct4dataset.zScorePre(rowCount,1) = CellmeanPre;
            struct4dataset.zScorePost(rowCount,1) = CellmeanPost;
            
            
            %Cellvariance_base=nanmean(pEvtMtx_base{i}{j}); % the variance is the squared stdev, so no need to square it later on.
            CellvariancePre=nanvar(Cellmean(:,data.preRange),[],2); %
            CellvariancePost=nanvar(Cellmean(:,data.postRange),[],2); %
            
            %sdPre=nanstd(CellvariancePre); %
            %sdCs=nanstd(CellvariancePost); %
            
            CellMeanDiff=abs(CellmeanPre-CellmeanPost); %
            CellSqrvariancePre=CellvariancePre; % the stddev^2 = the variance, no need to do anything
            CellSqrvariancePost=CellvariancePost; %
            CellSumSqrvariance=sqrt((CellSqrvariancePre+CellSqrvariancePost)/2); %
            
            %DPrime of Unit
            dSlctIdx{i}(1,j)=(CellMeanDiff)/(CellSumSqrvariance); % computes d: divides the absolute values of the mean differences by the sum of the squared deviations.
            %struct4dataset.DPrimeValue=dSlctIdx{i}(1,j);
            
            rowCount = rowCount+1;
        end;
    end;
    close(h);
    
   
    %% Create surrogate
    h = waitbar(0,'Initializing waitbar...');  % start waitbar to show progress
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Sur Mtx Init
    %tempSurMtx={size(pEvtMtx,2),1};
    %for i=1:sizepEvtMtx,2)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    for y=1:500; % this changes the number of times surrogates are run (use 1:1 if you dont need surrogate data).
        %     y
        PercentComplete = y/500;  waitbar(PercentComplete,h,sprintf('%f%% Complete...SURROGATE',PercentComplete));  %
        for i=1:size(pEvtMtx,2); %
            for j=1:size(pEvtMtx{i},2) %
                
                %% Scrambling data sets
                CellmeanOrig=zscore(nanmean(pEvtMtx{i}{j},1));
                x=CellmeanOrig(:);
                k=length(x);
                x=x(randperm(k));
                surMtx=reshape(x,size(CellmeanOrig));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % save the surrogate mtx for reproducability. 
                % need to add to futre versions
                %tempSurMtx{i}{j} = surMtx{i}{j}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                surMtx_pre=surMtx(:,data.preRange);
                surMtx_post=surMtx(:,data.postRange);
                
                
                CellmeanPre=nanmean(surMtx_pre,2);
                CellmeanPost=nanmean(surMtx_post,2);
                %                     Cellmean=nanmean(surMtx_post); %
                %                     CellmeanPre=nanmean(surMtx_pre(1:lLim)); %
                %                     CellmeanPost=nanmean(Cellmean(epoch)); %
                
                %Cellvariance=nanmean(surMtx_post); %
                %Cellvariance_pre=nanmean(surMtx_pre);
                CellvariancePre=nanvar(surMtx_pre,[],2); %
                CellvariancePost=nanvar(surMtx_post,[],2); %
                
                CellMeanDiff=abs(CellmeanPre-CellmeanPost); %
                CellSqrvariancePre=CellvariancePre; %
                CellSqrvariancePost=CellvariancePost; %
                CellSumSqrvariance=sqrt((CellSqrvariancePre+CellSqrvariancePost)/2); %
                surSlctIdx{y,i}(1,j)=(CellMeanDiff)/(CellSumSqrvariance); %
                
            end;
        end;
    end;
    close(h);
    
    %% plotting actual versus surrogate data - NEED TO REMOVE NANS AND INFS HERE; THEN NEED TO USE SURMTX TO GET STD AND MEAN TO GET CI
    
    surMtx=cell2mat(surSlctIdx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save the surrogate mtx for reproducability. 
    data.(curVar).DPrime.First3.surMtx = surMtx;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    surMtx_mean=nanmean(surMtx);
    surMtx_std=nanstd(surMtx);
    
    % GETTING THE UPPER CI FOR EACH NEURON AND THEN THE MEAN OF THE UPPER CI
    %CI=((surMtx_mean(:))+(surMtx_std(:)*2.575));%
    CI=((surMtx_mean(:))+(surMtx_std(:)*1.96));%
    CI_MEAN=nanmean(CI); % THIS WILL BE USED AS THRESHOLD FOR SELECTIVITY - CONSERVATIVE THRESHOLD
    
    % SORTING NEURONS FROM LOWEST DPRIME TO HIGHEST
    [sortData,Idx]=sort(cell2mat(dSlctIdx)); %
    
    %% NEED TO TAKE THE EACH NEURON THRESHOLD NOW INSTEAD OF ALL NEURON THRESHOLD AND COUNT THESE AS SIGNIFICNATLY SELECTIVE...
    CI_LOWER=((surMtx_mean(:))-(surMtx_std(:)*1.96));
    DP=cell2mat(dSlctIdx);
    DP=DP';
    CI_vs_DP=DP-CI;%% This is before all this insane sorting nonsense
    
    
    
    [sortCIvDP,CIDPIdx]=sort(CI_vs_DP);
    CellId=CellId(CIDPIdx); % SORT CELL ID BY SIGNIFICANCE
    
    % TRYING SOME PLOTS...
    %hold off;
    h=figure;
    plot(sortCIvDP,'ko'); hold on;
    plot(DP(CIDPIdx,:),'ro')
    x=zeros([1 60]); plot(x,'k:','LineWidth',2)
    plot(CI_LOWER(CIDPIdx,:),'b'); plot(CI(CIDPIdx,:),'b')
    xlim([0 length(sortData)]);
    plot(repmat(CI_MEAN,1,size(surMtx,2)),'g:','LineWidth',2);
    title 'Neuron Selectivity [sorted by significance]'
    
    % SIGNIFICANCE WITH INDIVIDUAL CI
    SIG=find(sortCIvDP>0);
    size(SIG)
    NOTSIG=find(sortCIvDP<=0);
    size(NOTSIG)
    
    (length(SIG))/(length(NOTSIG)+length(SIG))*(100);
    
    % SIGNIFICANCE WITH MEAN CI (OLD WAY)
    OLDSIG=find(DP>CI_MEAN);
    size(OLDSIG)
    OLDNOTSIG=find(DP<=CI_MEAN);
    size(OLDNOTSIG)
    
    (length(OLDSIG))/(length(OLDNOTSIG)+length(OLDSIG))*(100);
    
    %         rocMtx=cell2mat(iRoc);
    %         rocMtx=rocMtx(Idx);
    %         figure;
    %         scatter(sortData,rocMtx);
    
    %       cd 'D:\Users\Bailey\Documents\Dropbox\Mouse MEA\Mouse MEA\Mat\CL';
    %       save(sprintf('MkSel_Res_BAILEY__%02d',t));
    %        savefig(sprintf('FIG_NODRINK%d.fig',t),h); % will create FIG1, FIG2,...

    %clear all;
    %close all;
    
    %% using unsorted to get Direction etc
    struct4dataset.DPrimeSig(CI_vs_DP>0)=1;
    struct4dataset.DPrimeSig(~(CI_vs_DP>0))=0; %& ~isnan(CI_vs_DP)) = 0;
    struct4dataset.DPrimeCategory ...
        (struct4dataset.zScorePre-struct4dataset.zScorePost > 0 ...
        & struct4dataset.DPrimeSig) = {'decrease'};
    struct4dataset.DPrimeCategory ...
        (struct4dataset.zScorePre-struct4dataset.zScorePost < 0 ...
        & struct4dataset.DPrimeSig) = {'increase'};
        struct4dataset.DPrimeCategory ...
            (struct4dataset.DPrimeSig==0) = {'nochange'};
    
    data.(curVar).DPrime.First3.CellId = CellId;
    data.(curVar).DPrime.First3.SIG = SIG;
    data.(curVar).DPrime.First3.CIDPIdx = CIDPIdx;
    data.(curVar).DPrime.First3.dSlctIdx = dSlctIdx;
    
    data.(curVar).dataSet=struct2dataset(struct4dataset);
end


data.output.dataSet=[data.LL.dataSet; data.RL.dataSet];
data.output.dataSet.SessionType = nominal(data.output.dataSet.SessionType);
data.output.dataSet.Schedule = nominal(data.output.dataSet.Schedule);
data.output.dataSet.AnimalID = nominal(data.output.dataSet.AnimalID);
data.output.dataSet.DPrimeCategory = nominal(data.output.dataSet.DPrimeCategory);


%     struct4dataset.SessionType = cell(numRows,1);
%     struct4dataset.Schedule = cell(numRows,1);
%     struct4dataset.AnimalID = cell(numRows,1);
%     struct4dataset.CellID = cell(numRows,1);
%     struct4dataset.numEvents = nan(numRows,1);
%     struct4dataset.numAnalyzedEvents = nan(numRows,1);
%     struct4dataset.zScorePre = nan(numRows,1);
%     struct4dataset.zScorePost = nan(numRows,1);
%     struct4dataset.DPrimeSig = nan(numRows,1);
%     struct4dataset.DPrimeCategory = cell(numRows,1);
%     struct4dataset.zScoreOverall = nan(numRows,length(data.xA));

%% Here begins the old DPrimeDetails

%sort(data.LL04.DPrime.First3.CIDPIdx(data.LL04.DPrime.First3.SIG))
%% clean up the cells in CIDPIdx for LL and RL that do not have any firing during the period of interest
% VarList={'LL04','RL04','LL10','RL10'};
% for varIDX=[1 3]
%     toRemove=sort(unique( ...
%         [data.(VarList{varIDX}).DPrime.First3.CIDPIdx(max(data.(VarList{varIDX}).DPrime.First3.SIG)+1:end); ...
%         data.(VarList{varIDX+1}).DPrime.First3.CIDPIdx(max(data.(VarList{varIDX+1}).DPrime.First3.SIG)+1:end)]));
% end
%%

% may not need anything under here anymore
% 
% VarList={'LL','RL'};
% for varIDX=1:2    
%     tempZ=zscore(data.(VarList{varIDX}).MlMtx')';
%     SIG=data.(VarList{varIDX}).DPrime.First3.SIG;
%     CIDPIdx=data.(VarList{varIDX}).DPrime.First3.CIDPIdx;
%     newIDs=sort(CIDPIdx(SIG));
%     %     tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
%     %     [sorted, idx] = sort(mean( tempZ(:,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(:,data.xA>=0 & data.xA<=3),2 ));
%     diffmean= nanmean( tempZ(newIDs,data.preRange),2)- nanmean(tempZ(newIDs,data.postRange),2 );
%     data.(VarList{varIDX}).DPrime.First3.IncZ=tempZ(newIDs((diffmean<0)),:);
%     data.(VarList{varIDX}).DPrime.First3.DecZ=tempZ(newIDs((diffmean>0)),:);
%     data.(VarList{varIDX}).DPrime.First3.IncIDs=newIDs((diffmean<0));
%     data.(VarList{varIDX}).DPrime.First3.DecIDs=newIDs((diffmean>0));
%     data.(VarList{varIDX}).DPrime.First3.NCIDs=sort(CIDPIdx(~ismember(CIDPIdx,CIDPIdx(SIG))));
%     data.(VarList{varIDX}).DPrime.First3.NCZ=tempZ(data.(VarList{varIDX}).DPrime.First3.NCIDs);
%     diffmean=nanmean( tempZ(CIDPIdx,data.preRange),2)- nanmean(tempZ(CIDPIdx,data.postRange),2 );
%     diffdecidx=find(diffmean<0);
%     diffincidx=find(diffmean>0);
%     diffzerosidx=find(diffmean==0);
%     data.(VarList{varIDX}).DPrime.First3.overallSort=CIDPIdx([flipud(diffincidx);diffzerosidx;diffdecidx]);
% end
% 
% 
% %% New! Dprime Output for Stats
% VarList={'LL','RL'};
% %TimePointList={'E' 'E' 'L' 'L'};
% ScheduleList={'RI' 'RR'};
% %% Initialization
% data.Output.DPrime.Line.IDArrayInc=[data.LL.DPrime.First3.IncIDs;data.RL.DPrime.First3.IncIDs];
% data.Output.DPrime.Line.IDArrayDec=[data.LL.DPrime.First3.DecIDs;data.RL.DPrime.First3.DecIDs];
% data.Output.DPrime.Line.TimePointArrayInc=cell(length(data.Output.DPrime.Line.IDArrayInc),1);
% data.Output.DPrime.Line.TimePointArrayDec=cell(length(data.Output.DPrime.Line.IDArrayDec),1);
% data.Output.DPrime.Line.ScheduleArrayInc=cell(length(data.Output.DPrime.Line.IDArrayInc),1);
% data.Output.DPrime.Line.ScheduleArrayDec=cell(length(data.Output.DPrime.Line.IDArrayDec),1);
% data.Output.DPrime.Line.DataArrayInc=NaN(length(data.Output.DPrime.Line.IDArrayInc),length(data.xA));
% data.Output.DPrime.Line.DataArrayDec=NaN(length(data.Output.DPrime.Line.IDArrayDec),length(data.xA));
% length1=length(data.LL.DPrime.First3.CIDPIdx);
% length2=length(data.RL.DPrime.First3.CIDPIdx);
% data.Output.DPrime.Category.IDArray=[1:length1 length1+1:length1+length2]';
% %data.Output.DPrime.Category.TimePoint=cell(length1+length2,1);
% data.Output.DPrime.Category.Cat1=cell(length1+length2,1);
% data.Output.DPrime.Category.Cat2=cell(length1+length2,1);
% %data.Output.DPrime.Category.TimePoint(1:length1) = {'E'};
% %data.Output.DPrime.Category.TimePoint(length1+1:length1+length2) = {'L'};
% xInc=1;
% xDec=1;
% for varIDX=1:2
%     if varIDX<2
%         cellOffset=0;
%     else
%         cellOffset=length(data.LL.DPrime.First3.CIDPIdx);
%     end
%     if varIDX==2|| varIDX==4
%         CatIDX='Cat2';
%     else
%         CatIDX='Cat1';
%     end
%     data.Output.DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.IncIDs)={'I'};
%     data.Output.DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.DecIDs)={'D'};
%     data.Output.DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.NCIDs)={'N'};
%     
%     
%     tempInc=length(data.(VarList{varIDX}).DPrime.First3.IncIDs);
%     tempDec=length(data.(VarList{varIDX}).DPrime.First3.DecIDs);
%     % to delete
%     %data.Output.DPrime.tempInc=length(data.(VarList{varIDX}).DPrime.First3.IncIDs);
%     %data.Output.DPrime.tempDec=length(data.(VarList{varIDX}).DPrime.First3.DecIDs);
%     
%     %data.Output.DPrime.Line.TimePointArrayInc(xInc:xInc+tempInc-1)={TimePointList{varIDX}};
%     %data.Output.DPrime.Line.TimePointArrayDec(xDec:xDec+tempDec-1)={TimePointList{varIDX}};
%     data.Output.DPrime.Line.ScheduleArrayInc(xInc:xInc+tempInc-1)={ScheduleList{varIDX}};
%     data.Output.DPrime.Line.ScheduleArrayDec(xDec:xDec+tempDec-1)={ScheduleList{varIDX}};
%     data.Output.DPrime.Line.DataArrayInc(xInc:xInc+tempInc-1,:)=data.(VarList{varIDX}).DPrime.First3.IncZ;
%     data.Output.DPrime.Line.DataArrayDec(xDec:xDec+tempDec-1,:)=data.(VarList{varIDX}).DPrime.First3.DecZ;
%     xInc=xInc+tempInc;
%     xDec=xDec+tempDec;
% end
% %% get the zscore data for outputting for stats
% for varIDX=1:2
%     data.Output.zDataOutput{varIDX}=zscore(data.(VarList{varIDX}).MlMtx')';
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Getting the pie charts
% VarList={'LL04','RL04','LL10','RL10'};
% figure;
% for varIDX=1:4
%     subplot(2,2,varIDX);
%     pie(data.(VarList{varIDX}).DPrime.First3.pieData,{['Inc: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(1))] ...
%         ['NoChange: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(2))] ...
%         ['Dec: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(3))]});
% end
% %% getting the line graphs
% figure;
% CondList={'IncZ','DecZ'};
% meanRange=1:4;
% stderrRange=5:8;
% for condIDX=1:2
%     subplot(2,1,condIDX);
%     for varIDX = 1:4
%         y(varIDX,:)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
%         yStdErr(varIDX,:)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
%     end
%     plot(data.xA,y,'linewidth',2.5);
%     hold on
%     %plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','.','markersize',3,'linestyle','none');
%     plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
%     plot(downsample(data.xA,1),downsample(y+yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
%     %plot(data.xA,y+yStdErr,'linestyle',':');
%     title(['Overall zscore: ' CondList{condIDX}]);
%     axis([-3 3 -1.5 2.5]);
%     legend (VarList{:});
% end
%
%

% VarList={'LL04','RL04','LL10','RL10'};
% %PC = 4;
% figure;
% for varIDX=1:4
%     zData=zscore(data.(VarList{varIDX}).MlMtx')';
%     y(varIDX,:)=nanmean(zData);
%     yStdErr(varIDX,:)=nanstd(zData)/sqrt(size(zData,1));
% end
%
% %subplot(2,1,condIDX);
% plot(data.xA,y,'linewidth',2.5);
% hold on
% %plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','.','markersize',3,'linestyle','none');
% plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
% plot(downsample(data.xA,1),downsample(y+yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
% %plot(data.xA,y+yStdErr,'linestyle',':');
% title(['Overall zscore']);
% axis([-3 3 -1 1]);
% legend (VarList{:});
%
% %% DPrimeSorting Heatmap
% VarList={'LL04','RL04','LL10','RL10'};
% %PC = 4;
% figure;
% for varIDX=1:4
%     subplot(2,2,varIDX);
%     tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
%     imagesc(tempZ(data.(VarList{varIDX}).DPrime.First3.overallSort,:),[-2 2]);
%     title(['DPrime sorted: ' VarList{varIDX}]);
% end
% end
end


