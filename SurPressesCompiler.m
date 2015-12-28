VarList={'LL04','RL04','LL10','RL10'};
temp=[];
totCount=1;
%dataSet='SurSpike';
for varIDX=1:4
    for CellIDX=1:length(data.(VarList{varIDX}).SurSpike)
        multVar=length(data.(VarList{varIDX}).pEvt{CellIDX});%%%no pevt!!!!!!
        numTrials=length(data.(VarList{varIDX}).SurSpike{CellIDX});
        curRange=totCount:totCount+numTrials-1;
        temp(curRange,1)=varIDX;%CondID
        temp(curRange,2)=CellIDX;%AnimalID
        temp(curRange,3)=multVar;%CellCount
        temp(curRange,4)=curRange-totCount+1;%TrialID
        %Rf2BB
        temp(curRange,5:6)= [data.(VarList{varIDX}).SurSpike{CellIDX}];
        %0=Bad 1=Good
        temp(curRange,7)=~isnan(data.(VarList{varIDX}).pEvt{CellIDX}{1}(:,1));
        totCount = totCount + numTrials;
    end
end
for x=1:4
    tempmean=mean(temp(temp(:,1)==x,5:6),1);
    temperr=std(temp(temp(:,1)==x,5:6),1)/sqrt(length(temp(temp(:,1)==x ,5)));
    meanVals(x,[1 3])=tempmean;
    meanVals(x,[2 4])=temperr;
    %    meanValsString{x}= [num2str(tempmean(1),'%.2f') char(177) num2str(temperr,'%.2f')];
    meanValsUsed(x,[1 3])=mean(temp(temp(:,1)==x & temp(:,end)==1,5:6),1);
    meanValsUsed(x,[2 4])=std(temp(temp(:,1)==x & temp(:,end)==1 ,5:6),1)/sqrt(length(temp(temp(:,1)==x & temp(:,end)==1 ,5)));
end
meanVals=reshape(meanVals',1,16);
for x=1:8
    meanValsString{x}= [num2str(meanVals((x*2)-1),'%.2f') setstr(177) num2str(meanVals(x*2),'%.2f')];
end
meanValsUsed=reshape(meanValsUsed',1,16);
for x=1:8
    meanValsUsedString{x}= [num2str(meanValsUsed((x*2)-1),'%.2f') setstr(177) num2str(meanValsUsed(x*2),'%.2f')];
end
for x=1:4
    tempanimalMean=[];
    for y=1:max(temp(temp(:,1)==x,2))
        %meanVals(x,:)=mean(temp(temp(:,1)==x,5:6),1);
        %meanVals(x+4,:)=std(temp(temp(:,1)==x,5:6),1)/sqrt(length(temp(temp(:,1)==x ,5)));
        tempanimalMean(y,:)=mean(temp(temp(:,1)==x & temp(:,end)==1 & temp(:,2)==y,5:6),1);
        %meanValsUsed(x+4,:)=std(temp(temp(:,1)==x & temp(:,end)==1 ,5:6),1)/sqrt(length(temp(temp(:,1)==x & temp(:,end)==1 ,5)));
    end
    meanValsByAnimal(x,[1 3])=mean(tempanimalMean,1);
    meanValsByAnimal(x,[2 4])=std(tempanimalMean,1)/sqrt(length(tempanimalMean));
end
meanValsByAnimal=reshape(meanValsByAnimal',1,16);
for x=1:8
    meanValsByAnimalString{x}= [num2str(meanValsByAnimal((x*2)-1),'%.2f') setstr(177) num2str(meanValsByAnimal(x*2),'%.2f')];
end
%% Get nuber of cells in groups
tempData=nan(8,6);

VarList={'LL04','RL04','LL10','RL10'};
for varIDX=1:4
    tempCount = 1;

    cellCounts = cellfun(@length,data.(VarList{varIDX}).pEvt);
    for x=1:length(data.(VarList{varIDX}).pEvt)
        tempData(varIDX,x)=length(intersect (data.(VarList{varIDX}).DPrime.First3.IncIDs,tempCount:tempCount+cellCounts(x)-1));
        tempData(varIDX+4,x)=length(intersect (data.(VarList{varIDX}).DPrime.First3.DecIDs,tempCount:tempCount+cellCounts(x)-1));
        tempCount=tempCount+cellCounts(x)-1;
    end
end