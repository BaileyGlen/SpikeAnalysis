

VarList={'LL04','RL04','LL10','RL10'};
temp=[ 0    0.4470    0.7410;...
0.8500    0.3250    0.0980;...
0.9290    0.6940    0.1250;...
0.4940    0.1840    0.5560];
newcmap = [temp(1,:);temp(3,:);temp(4,:);temp(2,:)];

%% overal means zscored - plotted as 1 plot

VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;

for varIDX=1:4
    zData=zscore(data.(VarList{varIDX}).MlMtx')';
    y(:,varIDX)=nanmean(zData);
    yStdErr(:,1,varIDX)=nanstd(zData)/sqrt(size(zData,1));
end
%%

fighandle=figure;
%newStdErr=permute(yStdErr,[1,2,3]);
[hl, hp]=boundedline(data.xA, y, yStdErr ,'cmap',newcmap, 'alpha');
%title(['Overall zscore']);
axis([-3 3 -.6 .8]);
%legend (VarList{:});
set(hl,'linewidth',2.5);
set(fighandle, 'Position', [100, 100, 350, 250]);


%     %subplot(2,1,condIDX);
%     plot(data.xA,y,'linewidth',2.5);
%     hold on
%     %plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','.','markersize',3,'linestyle','none');
%     plot(downsample(data.xA,1),downsample(y-yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
%     plot(downsample(data.xA,1),downsample(y+yStdErr,1),'marker','none','markersize',3,'linestyle',':','LineWidth',.5);
%     %plot(data.xA,y+yStdErr,'linestyle',':');
%     title(['Overall zscore']);
%     axis([-3 3 -1 1]);
%     legend (VarList{:});
%% Inc and Dec
fighandle=figure;
CondList={'IncZ','DecZ'};
meanRange=1:4;
stderrRange=5:8;
for condIDX=1:2
    subplot(2,1,condIDX);
    for varIDX = 1:4
        y(:,varIDX)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
        yStdErr(:,1,varIDX)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
    end
    [hl, hp]=boundedline(data.xA, y, yStdErr ,'cmap',newcmap, 'alpha');
    set(fighandle, 'Position', [100, 100, 350, 250]);
    if condIDX==1
        tmpStr='Increasing';
    else
        tmpStr='Decreasing';
    end
    title(['Reinforcer Consumption Modulated Cells: ' tmpStr]);
    set(hl,'linewidth',2.5);
    plot([-3 13],[0 0],'linestyle','--','color',[.5 .5 .5],'linewidth',1);
    axis([-3 5 -1.5 2.5]);
    
    %h_legend=legend (VarList{:});
    %set(h_legend,'FontSize',8,'Orientation','horizontal','Location','best') 
    
end

%%heatmap
%% DPrimeSorting Heatmap
VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;
fighandle=figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    set(fighandle, 'Position', [100, 100, 350, 250]);
    imagesc(data.xA,1:size(tempZ,1),tempZ(data.(VarList{varIDX}).DPrime.First3.CIDPIdx(data.(VarList{varIDX}).DPrime.First3.overallSort),:),[-2 2]);
    %title(['DPrime sorted: ' VarList{varIDX}]);
    axis([-3 5 1 size(tempZ,1)]);
    set(gca,'YTick',[1 size(tempZ,1)])
end

%% Getting the pie charts
VarList={'LL04','RL04','LL10','RL10'};
fighandle=figure;
for varIDX=1:4
    subplot(2,2,varIDX);
    set(fighandle, 'Position', [100, 100, 350, 250]);
    pie(data.(VarList{varIDX}).DPrime.First3.pieData);
%     ,{['Inc: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(1))] ...
%         ['NoChange: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(2))] ...
%         ['Dec: ' num2str(data.(VarList{varIDX}).DPrime.First3.pieData(3))]});
end