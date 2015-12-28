VarList={'LL04','RL04','LL10','RL10'};
temp=[];
totCount=1;
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