function [ data ] = DPrimeDetails( data ,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if type=='C' %consumption
    %these are reversed pre vs post in order to allow inc and dec to stay
    %the same
    preRange=data.xA>=8 & data.xA<=13;
    postRange=data.xA>=0 & data.xA<=5;
elseif type=='P'
    preRange=data.xA>=-3 & data.xA<=0;
    postRange=data.xA>=0 & data.xA<=3;
else error('myApp:argChk','Invalid Input for "Type"')
end
%sort(data.LL04.DPrime.First3.CIDPIdx(data.LL04.DPrime.First3.SIG))
%% clean up the cells in CIDPIdx for LL and RL that do not have any firing during the period of interest
% VarList={'LL04','RL04','LL10','RL10'};
% for varIDX=[1 3]
%     toRemove=sort(unique( ...
%         [data.(VarList{varIDX}).DPrime.First3.CIDPIdx(max(data.(VarList{varIDX}).DPrime.First3.SIG)+1:end); ...
%         data.(VarList{varIDX+1}).DPrime.First3.CIDPIdx(max(data.(VarList{varIDX+1}).DPrime.First3.SIG)+1:end)]));
% end
%%
VarList={'LL04','RL04','LL10','RL10'};
for varIDX=1:4
    tempZ=zscore(data.(VarList{varIDX}).MlMtx')';
    SIG=data.(VarList{varIDX}).DPrime.First3.SIG;
    CIDPIdx=data.(VarList{varIDX}).DPrime.First3.CIDPIdx;
    newIDs=sort(CIDPIdx(SIG));
    %     tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    %     [sorted, idx] = sort(mean( tempZ(:,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(:,data.xA>=0 & data.xA<=3),2 ));
    diffmean=mean( tempZ(newIDs,preRange),2)- mean(tempZ(newIDs,postRange),2 );
    data.(VarList{varIDX}).DPrime.First3.IncZ=tempZ(newIDs((diffmean<0)),:);
    data.(VarList{varIDX}).DPrime.First3.DecZ=tempZ(newIDs((diffmean>0)),:);
    data.(VarList{varIDX}).DPrime.First3.IncIDs=newIDs((diffmean<0));
    data.(VarList{varIDX}).DPrime.First3.DecIDs=newIDs((diffmean>0));
    data.(VarList{varIDX}).DPrime.First3.NCIDs=sort(CIDPIdx(~ismember(CIDPIdx,CIDPIdx(SIG))));
    data.(VarList{varIDX}).DPrime.First3.NCZ=tempZ(data.(VarList{varIDX}).DPrime.First3.NCIDs);
    diffmean=mean( tempZ(CIDPIdx,preRange),2)- mean(tempZ(CIDPIdx,postRange),2 );
    diffdecidx=find(diffmean<0);
    diffincidx=find(diffmean>0);
    diffzerosidx=find(diffmean==0);
    data.(VarList{varIDX}).DPrime.First3.overallSort=CIDPIdx([flipud(diffincidx);diffzerosidx;diffdecidx]);
end
%% New! Dprime Output for Stats
VarList={'LL04','RL04','LL10','RL10'};
TimePointList={'E' 'E' 'L' 'L'};
ScheduleList={'RI' 'RR' 'RI' 'RR'};
data.Output.DPrime.Line.IDArrayInc=[data.LL04.DPrime.First3.IncIDs;data.RL04.DPrime.First3.IncIDs;...
    data.LL10.DPrime.First3.IncIDs+51;data.RL10.DPrime.First3.IncIDs+51];
data.Output.DPrime.Line.IDArrayDec=[data.LL04.DPrime.First3.DecIDs;data.RL04.DPrime.First3.DecIDs;...
    data.LL10.DPrime.First3.DecIDs+51;data.RL10.DPrime.First3.DecIDs+51];
data.Output.DPrime.Line.TimePointArrayInc=cell(length(data.Output.DPrime.Line.IDArrayInc),1);
data.Output.DPrime.Line.TimePointArrayDec=cell(length(data.Output.DPrime.Line.IDArrayDec),1);
data.Output.DPrime.Line.ScheduleArrayInc=cell(length(data.Output.DPrime.Line.IDArrayInc),1);
data.Output.DPrime.Line.ScheduleArrayDec=cell(length(data.Output.DPrime.Line.IDArrayDec),1);
data.Output.DPrime.Line.DataArrayInc=NaN(length(data.Output.DPrime.Line.IDArrayInc),length(data.xA));
data.Output.DPrime.Line.DataArrayDec=NaN(length(data.Output.DPrime.Line.IDArrayDec),length(data.xA));
length1=length(data.LL04.DPrime.First3.CIDPIdx);
length2=length(data.LL10.DPrime.First3.CIDPIdx);
data.Output.DPrime.Category.IDArray=[1:length1 length1+1:length1+length2]';
data.Output.DPrime.Category.TimePoint=cell(length1+length2,1);
data.Output.DPrime.Category.Cat1=cell(length1+length2,1);
data.Output.DPrime.Category.Cat2=cell(length1+length2,1);
data.Output.DPrime.Category.TimePoint(1:length1) = {'E'};
data.Output.DPrime.Category.TimePoint(length1+1:length1+length2) = {'L'};
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
    data.Output.DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.IncIDs+cellOffset)={'I'};
    data.Output.DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.DecIDs+cellOffset)={'D'};
    data.Output.DPrime.Category.(CatIDX)(data.(VarList{varIDX}).DPrime.First3.NCIDs+cellOffset)={'N'};
    
    
    tempInc=length(data.(VarList{varIDX}).DPrime.First3.IncIDs);
    tempDec=length(data.(VarList{varIDX}).DPrime.First3.DecIDs);
    % to delete
    %data.Output.DPrime.tempInc=length(data.(VarList{varIDX}).DPrime.First3.IncIDs);
    %data.Output.DPrime.tempDec=length(data.(VarList{varIDX}).DPrime.First3.DecIDs);
    
    data.Output.DPrime.Line.TimePointArrayInc(xInc:xInc+tempInc-1)={TimePointList{varIDX}};
    data.Output.DPrime.Line.TimePointArrayDec(xDec:xDec+tempDec-1)={TimePointList{varIDX}};
    data.Output.DPrime.Line.ScheduleArrayInc(xInc:xInc+tempInc-1)={ScheduleList{varIDX}};
    data.Output.DPrime.Line.ScheduleArrayDec(xDec:xDec+tempDec-1)={ScheduleList{varIDX}};
    data.Output.DPrime.Line.DataArrayInc(xInc:xInc+tempInc-1,:)=data.(VarList{varIDX}).DPrime.First3.IncZ;
    data.Output.DPrime.Line.DataArrayDec(xDec:xDec+tempDec-1,:)=data.(VarList{varIDX}).DPrime.First3.DecZ;
    xInc=xInc+tempInc;
    xDec=xDec+tempDec;
end
%% get the zscore data for outputting for stats
for varIDX=1:4
    data.Output.zDataOutput{varIDX}=zscore(data.(VarList{varIDX}).MlMtx')';
end
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