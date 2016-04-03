

VarList={'LL','RL'};
temp=[ 0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560];
newcmap = [temp(1,:);temp(3,:);temp(4,:);temp(2,:)];

%% overal means zscored - plotted as 1 plot

%VarList={'LL04','RL04','LL10','RL10'};
%PC = 4;

for varIDX=1:2
    zData=zscore(data.(VarList{varIDX}).MlMtx')';
    y(:,varIDX)=nanmean(zData);
    yStdErr(:,1,varIDX)=nanstd(zData)/sqrt(size(zData,1));
end

fighandle=figure;
%newStdErr=permute(yStdErr,[1,2,3]);
[hl, hp]=boundedline(data.xA, y, yStdErr ,'cmap',newcmap, 'alpha');
%title(['Overall zscore']);
axis([-2.8 8 -.75 .95]);
%legend (VarList{:});
set(hl,'linewidth',2.5);
set(fighandle, 'Position', [100, 100, 420, 208]);
set(gca,'YTick', [-.4 0 .4 .8]);
set(gca,'XTick', [-2 0 2 4 6 8]);
set(gca,'FontSize',8)
set(gca,'FontName','Calibri')
set(gca,'Linewidth',.75);
axis on;
box off;


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
%export_fig Fig2A.jpg -m3.5 -transparent
%export_fig (fighandle, 'tempTitle.jpg');
%% Inc and Dec
% Current Working!!!!

CondList={'IncZ','DecZ'};
%TimeList={'Early','Extended'};
y=[];
yStdErr=[];
%newcmap = [temp(4,:);temp(2,:)];
newcmap = [temp(1,:);temp(3,:);temp(4,:);temp(2,:)];
for condIDX = 1:2
    
    
    %set(fighandle, 'Position', [100, 100, 202, 215]);
    
    for varIDX = 1:2
        y(:,varIDX)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
        yStdErr(:,1,varIDX)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
        %y(:,varIDX-2)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
        %yStdErr(:,1,varIDX-2)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
    end
    
    fighandle=figure ('position',[100, 100, 200, 200]);%, 'outerposition',[90, 90, 250, 200]);
    %curAxes=axes('position', [0.1300    0.1100+(.5*x-1)    0.7750    0.3412]);
    %axisPos=get(gca, 'position');
    %set (gca,'position', [0.1300+(    0.1100    0.7750    0.3412]);
    %tempCmap=newcmap(timeIDX*2-1:timeIDX*2,:);
    subplot('Position',[.2 .2 .63 .59]);
    
    [hl, hp]=boundedline(data.xA, y(), yStdErr() ,'cmap',newcmap);
    %set(gca,'position', [0.1300    0.1100+(.5*x-1)    0.7750    0.3412]);
    %axisPos=get(gca, 'position');
    %set(gca, 'position', [axisPos(1)-(~mod(condIDX,2)*.04) axisPos(2)+.07 axisPos(3) axisPos(4)]);
    
    %     if condIDX==1
    %         tmpStr='Increasing';
    %     else
    %         tmpStr='Decreasing';
    %     end
    %title(['Reinforcer Consumption Modulated Cells: ' tmpStr]);
    set(hl,'linewidth',2.2);
    plot([-3 13],[0 0],'linestyle','-','color',[.5 .5 .5],'linewidth',.9);

    axis([-2.8 8 -1.5 2.5]);
    set(gca,'YTick', [-1 0 1 2]);
    set(gca,'XTick', [-2 0 2]);
    set(gca,'FontSize',8)
    set(gca,'FontName','Calibri')
    set(gca,'Linewidth',.75);
    
    
    %h_legend=legend (VarList{:});
    %set(h_legend,'FontSize',8,'Orientation','horizontal','Location','best')
end    
    



%% DPrimeSorting Heatmap
VarList={'LL','RL'};
%PC = 4;
fighandle=figure;
for varIDX=1:2
    set(gca,'DefaultTextFontSize',18+varIDX*2)
    subplot(1,2,varIDX);
    tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    set(fighandle, 'Position', [100, 100, 280, 208]);
    %set(gca,'DefaultTextFontSize',18+varIDX*2)
    imagesc(data.xA,1:size(tempZ,1),tempZ(data.(VarList{varIDX}).DPrime.First3.overallSort,:),[-2 2]);
    %title(['DPrime sorted: ' VarList{varIDX}]);
    axis([-3 3 1 size(tempZ,1)]);
    axisPos=get(gca, 'position');
    set(gca, 'position', [axisPos(1)-(~mod(varIDX,2)*.02) axisPos(2)+((varIDX<=2)*.02) axisPos(3)/1.08 axisPos(4)]);
    set(gca,'YTick',[1 size(tempZ,1)])
    set(gca,'FontSize',8)
    set(gca,'FontName','Calibri')
end
axes('cLim',[-2 2],'Position', [.1 .1 .9 .9], 'Visible', 'off');
cb=colorbar('location', 'East','YAxisLocation','right');
%set(cbar_handle, 'YAxisLocation','right')
cbPos=get(cb, 'position');
set(cb, 'position', [cbPos(1)+.03 cbPos(2)-.025 cbPos(3)/5 cbPos(4)]);
set(cb,'FontSize',8)
set(cb,'FontName','Calibri')
export_fig RfDelHeat.png -m3.5 -transparent
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


%% Inc and Dec --- UnRf
% Current Working!!!!

CondList={'IncZ','DecZ'};
TimeList={'Early','Extended'};
meanRange=1:2;
stderrRange=3:4;
y=[];
yStdErr=[];
%newcmap = [temp(4,:);temp(2,:)];
newcmap = [temp(1,:);temp(3,:);temp(4,:);temp(2,:)];
for condIDX=1:2
    

    %set(fighandle, 'Position', [100, 100, 202, 215]);
    
    for varIDX = 1:2
        y(:,varIDX)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
        yStdErr(:,1,varIDX)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
        %y(:,varIDX-2)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
        %yStdErr(:,1,varIDX-2)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
    end
    for timeIDX=1:1
        fighandle=figure ('position',[100, 100, 250, 125]);%, 'outerposition',[90, 90, 250, 200]);
        %curAxes=axes('position', [0.1300    0.1100+(.5*x-1)    0.7750    0.3412]);
        %axisPos=get(gca, 'position');
        %set (gca,'position', [0.1300+(    0.1100    0.7750    0.3412]);
        tempCmap=newcmap(timeIDX*2-1:timeIDX*2,:);
        subplot('Position',[.2 .2 .63 .59]);
        [hl, hp]=boundedline(data.xA, y(:,timeIDX*2-1:timeIDX*2), yStdErr(:,timeIDX*2-1:timeIDX*2) ,'cmap',tempCmap);
        %set(gca,'position', [0.1300    0.1100+(.5*x-1)    0.7750    0.3412]);
        %axisPos=get(gca, 'position');
        %set(gca, 'position', [axisPos(1)-(~mod(condIDX,2)*.04) axisPos(2)+.07 axisPos(3) axisPos(4)]);
        
        %     if condIDX==1
        %         tmpStr='Increasing';
        %     else
        %         tmpStr='Decreasing';
        %     end
        %title(['Reinforcer Consumption Modulated Cells: ' tmpStr]);
        set(hl,'linewidth',2.2);
        plot([-3 13],[0 0],'linestyle','-','color',[.5 .5 .5],'linewidth',.9);
        %%This is for adding sig lines
%         if condIDX==2
%             if timeIDX==1 plot([0.4 0.7],[1 1],'linestyle','-','color',[0 0 0],'linewidth',2);
%             else plot([0.9 1.3],[.2 .2],'linestyle','-','color',[0 0 0],'linewidth',2);
%             end
%         end
        axis([-2.8 2.8 -1.5 2.5]);
        set(gca,'YTick', [-1 0 1 2]);
        set(gca,'XTick', [-2 0 2 4 6]);
        set(gca,'FontSize',8)
        set(gca,'FontName','Calibri')
        set(gca,'Linewidth',.75);    
        if condIDX==1 tempTitle=['UnRfInc' TimeList{timeIDX} '.png'];
        else tempTitle=['UnRfDec' TimeList{timeIDX} '.png'];
        end
        
        export_fig (fighandle, '-m3.5', '-transparent',tempTitle);
        
        
        %h_legend=legend (VarList{:});
        %set(h_legend,'FontSize',8,'Orientation','horizontal','Location','best')
    end


end

%% Inc and Dec -- UnRf Simple Inc Dec
fighandle=figure;
CondList={'IncZ','DecZ'};
meanRange=1:4;
stderrRange=5:8;
y=[];
yStdErr=[];
%newcmap = [temp(4,:);temp(2,:)];
newcmap = [temp(1,:);temp(3,:);temp(4,:);temp(2,:)];
set(fighandle, 'Position', [100, 100, 350, 200]);
for condIDX=1:2
    subplot(2,1,condIDX);
    for varIDX = 1:4
        y(:,varIDX)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
        yStdErr(:,1,varIDX)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
        %y(:,varIDX-2)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
        %yStdErr(:,1,varIDX-2)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
        axisPos=get(gca, 'position');
        
    end
    [hl, hp]=boundedline(data.xA, y, yStdErr ,'cmap',newcmap, 'alpha');
    %axisPos=get(gca, 'position');
    %set(gca, 'position', [axisPos(1)-(~mod(condIDX,2)*.04) axisPos(2)+.07 axisPos(3) axisPos(4)]);
    
    if condIDX==1
        tmpStr='Increasing';
    else
        tmpStr='Decreasing';
    end
    %title(['Reinforcer Consumption Modulated Cells: ' tmpStr]);
    set(hl,'linewidth',2.2);
    plot([-3 13],[0 0],'linestyle','--','color',[.5 .5 .5],'linewidth',.9);
    axis([-2.8 2.8 -1.5 2.5]);
    set(gca,'YTick', [-1 0 1 2]);
    set(gca,'XTick', [-2 0 2 4 6]);
    set(gca,'FontSize',8)
    set(gca,'FontName','Calibri')
    set(gca,'Linewidth',.75);
    
    %h_legend=legend (VarList{:});
    %set(h_legend,'FontSize',8,'Orientation','horizontal','Location','best')
    
end
export_fig UnRfD.png -m3.5 -transparent;
%% Inc and Dec -- RfDel Early vs Extended 

CondList={'IncZ','DecZ'};
meanRange=1:2;
stderrRange=3:4;

%newcmap = [temp(4,:);temp(2,:)];
newcmap = [temp(1,:);temp(3,:);temp(4,:);temp(2,:)];

for condIDX=1:2
    fighandle=figure;
    set(fighandle, 'Position', [100, 100, 250, 210]);
    %subplot(2,1,condIDX);
    %fighandle=figure;
    y=[];
    yStdErr=[];
    for varIDX=(1:2)
        zData=zscore(data.(VarList{varIDX+(2*(condIDX-1))}).MlMtx')';
        y(:,varIDX)=nanmean(zData);
        yStdErr(:,1,varIDX)=nanstd(zData)/sqrt(size(zData,1));
        axisPos=get(gca, 'position');
    end
    
    %     for varIDX = 1:4
%         y(:,varIDX)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
%         yStdErr(:,1,varIDX)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
%         %y(:,varIDX-2)=mean(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1);
%         %yStdErr(:,1,varIDX-2)=std(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1)/sqrt(size(data.(VarList{varIDX}).DPrime.First3.(CondList{condIDX}),1));
%         
%         
%     end
    [hl, hp]=boundedline(data.xA, y, yStdErr ,'cmap',newcmap((1:2)+(2*(condIDX-1)),:), 'alpha');
    %axisPos=get(gca, 'position');
    %set(gca, 'position', [axisPos(1)-(~mod(condIDX,2)*.04) axisPos(2)+.07 axisPos(3) axisPos(4)]);
    
    if condIDX==1
        tmpStr='Increasing';
    else
        tmpStr='Decreasing';
    end
    %title(['Reinforcer Consumption Modulated Cells: ' tmpStr]);
    set(hl,'linewidth',2.2);
    plot([-3 13],[0 0],'linestyle','--','color',[.5 .5 .5],'linewidth',.9);
    axis([-2 4 -.5 .8]);
    set(gca,'YTick', [-.4 0 .4 .8]);
    set(gca,'XTick', [-2 0 2 4 6]);
    set(gca,'FontSize',8)
    set(gca,'FontName','Calibri')
    set(gca,'Linewidth',.75);
    
    %h_legend=legend (VarList{:});
    %set(h_legend,'FontSize',8,'Orientation','horizontal','Location','best')
    figTitle=['RfDel_' num2str(condIDX) '.png'];
    export_fig (fighandle, '-m3.5', '-transparent',figTitle);
end
