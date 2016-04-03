VarList={'LL','RL'};
data.data.preRange=data.xA>=-3 & data.xA<=0;
data.data.postRange=data.xA>=0 & data.xA<=3;
for varIDX=1:2    
    tempZ=zscore(data.(VarList{varIDX}).MlMtx')';
    SIG=data.(VarList{varIDX}).DPrime.First3.SIG;
    CIDPIdx=data.(VarList{varIDX}).DPrime.First3.CIDPIdx;
    newIDs=sort(CIDPIdx(SIG));
    %     tempZ= zscore(data.(VarList{varIDX}).MlMtx')';
    %     [sorted, idx] = sort(mean( tempZ(:,data.xA>=-3 & data.xA<=0),2)- mean(tempZ(:,data.xA>=0 & data.xA<=3),2 ));
    diffmean= nanmean( tempZ(newIDs,data.preRange),2)- nanmean(tempZ(newIDs,data.postRange),2 );
    data.(VarList{varIDX}).DPrime.First3.IncZ=tempZ(newIDs((diffmean<0)),:);
    data.(VarList{varIDX}).DPrime.First3.DecZ=tempZ(newIDs((diffmean>0)),:);
    data.(VarList{varIDX}).DPrime.First3.IncIDs=newIDs((diffmean<0));
    data.(VarList{varIDX}).DPrime.First3.DecIDs=newIDs((diffmean>0));
    data.(VarList{varIDX}).DPrime.First3.NCIDs=sort(CIDPIdx(~ismember(CIDPIdx,CIDPIdx(SIG))));
    data.(VarList{varIDX}).DPrime.First3.NCZ=tempZ(data.(VarList{varIDX}).DPrime.First3.NCIDs);
    diffmean=nanmean( tempZ(CIDPIdx,data.preRange),2)- nanmean(tempZ(CIDPIdx,data.postRange),2 );
    diffdecidx=find(diffmean<0);
    diffincidx=find(diffmean>0);
    diffzerosidx=find(diffmean==0);
    data.(VarList{varIDX}).DPrime.First3.overallSort=CIDPIdx([flipud(diffincidx);diffzerosidx;diffdecidx]);
end