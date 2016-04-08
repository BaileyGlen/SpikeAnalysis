function [PCA] = ExaminePCA(  )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Run The PCA
mnSpk = 0;
mxSpk=1000;
rasterBin=.1;      % Size of bines for raster   
preEvt=3;            % time prior to Event in sec
postEvt=10;           % time post Event in sec
PLOT=0;
nmPC=4;
tempRange = [1:((preEvt+postEvt)/rasterBin)+1];
PCA.xA=[-1*preEvt:rasterBin:postEvt];
Day=4;
bEvt1={'LL_R'};

[Mtx, MlMtx, BlMtx, pEvt, pEvt_base]=mkPeriEvt_04(bEvt1,'r',Day,'zNo',PLOT, preEvt, postEvt, rasterBin);
xA=[-1*preEvt:rasterBin:postEvt];
k=find(sum(MlMtx')>mnSpk & sum(MlMtx')<mxSpk); % get rid of low FR ST
[coeff,score,latent,tsquared,explained,mu]=pca(bsxfun(@rdivide, MlMtx(k,tempRange),nanmean(BlMtx(k,:),2)));
%[coeff,score,latent,tsquared,explained,mu]=pca(MlMtx(k,:));
PCA.LL04.coeff=coeff;
PCA.LL04.score=score;
PCA.LL04.latent=latent;
PCA.LL04.tsquared=tsquared;
PCA.LL04.explained=explained;
PCA.LL04.mu=mu;
PCA.LL04.MlMtx = MlMtx;
PCA.LL04.BlMtx = BlMtx;
PCA.LL04.pEvt = pEvt;
PCA.LL04.pEvt_base = pEvt_base;

Day=4;
bEvt1={'RL_R'};

[Mtx, MlMtx, BlMtx, pEvt, pEvt_base]=mkPeriEvt_04(bEvt1,'r',Day,'zNo',PLOT, preEvt, postEvt, rasterBin);
xA=[-1*preEvt:rasterBin:postEvt];
k=find(sum(MlMtx')>mnSpk & sum(MlMtx')<mxSpk); % get rid of low FR ST
%[coeff,score,latent,tsquared,explained,mu]=pca(MlMtx(k,:));
[coeff,score,latent,tsquared,explained,mu]=pca(bsxfun(@rdivide, MlMtx(k,tempRange),nanmean(BlMtx(k,:),2)));
PCA.RL04.coeff=coeff;
PCA.RL04.score=score;
PCA.RL04.latent=latent;
PCA.RL04.tsquared=tsquared;
PCA.RL04.explained=explained;
PCA.RL04.mu=mu;
PCA.RL04.MlMtx = MlMtx;
PCA.RL04.BlMtx = BlMtx;
PCA.RL04.pEvt = pEvt;
PCA.RL04.pEvt_base = pEvt_base;

Day=10;
bEvt1={'LL_R'};

[Mtx, MlMtx, BlMtx, pEvt, pEvt_base]=mkPeriEvt_04(bEvt1,'r',Day,'zNo',PLOT, preEvt, postEvt, rasterBin);
xA=[-1*preEvt:rasterBin:postEvt];
k=find(sum(MlMtx')>mnSpk & sum(MlMtx')<mxSpk); % get rid of low FR ST
%[coeff,score,latent,tsquared,explained,mu]=pca(MlMtx(k,:));
[coeff,score,latent,tsquared,explained,mu]=pca(bsxfun(@rdivide, MlMtx(k,tempRange),nanmean(BlMtx(k,:),2)));
PCA.LL10.coeff=coeff;
PCA.LL10.score=score;
PCA.LL10.latent=latent;
PCA.LL10.tsquared=tsquared;
PCA.LL10.explained=explained;
PCA.LL10.mu=mu;
PCA.LL10.MlMtx = MlMtx;
PCA.LL10.BlMtx = BlMtx;
PCA.LL10.pEvt = pEvt;
PCA.LL10.pEvt_base = pEvt_base;

Day=10;
bEvt1={'RL_R'};

[Mtx, MlMtx, BlMtx, pEvt, pEvt_base]=mkPeriEvt_04(bEvt1,'r',Day,'zNo',PLOT, preEvt, postEvt, rasterBin);
xA=[-1*preEvt:rasterBin:postEvt];
k=find(sum(MlMtx')>mnSpk & sum(MlMtx')<mxSpk); % get rid of low FR ST
%[coeff,score,latent,tsquared,explained,mu]=pca(MlMtx(k,:));
[coeff,score,latent,tsquared,explained,mu]=pca(bsxfun(@rdivide, MlMtx(k,tempRange),nanmean(BlMtx(k,:),2)));
PCA.RL10.coeff=coeff;
PCA.RL10.score=score;
PCA.RL10.latent=latent;
PCA.RL10.tsquared=tsquared;
PCA.RL10.explained=explained;
PCA.RL10.mu=mu;
PCA.RL10.MlMtx = MlMtx;
PCA.RL10.BlMtx = BlMtx;
PCA.RL10.pEvt = pEvt;
PCA.RL10.pEvt_base = pEvt_base;
end

% figure;
% maxcond=4;
% if ~isfield(sortStruct,'conIDX')
%     sortStruct.conIDX=1;
%     sortStruct.data=cell(4);
%     sortStruct.score=cell(4);
%     sortStruct.sortIDX=cell(4);
%     display('Created the sortStruct');
% else
%     if sortStruct.conIDX>maxcond;
%         sortStruct = struct();
%         sortStruct.idx=1;
%         sortStruct.data=cell(4);
%         sortStruct.score=cell(4);
%         sortStruct.sortIDX=cell(4);
%         display('Reset the sortStruct');
%     end
% end
% for x=1:4
%     [PCA.score{x,PCA.conIDX},PCA.sortIDX{x,PCA.conIDX}]=SortPCA(score,x);
%     PCA.data{x,PCA.conIDX}=data(PCA.sortIDX{x,PCA.conIDX},:);
%     subplot(2,2,x);
%     imagesc(xA,1:size(PCA.data,1),zscore(PCA.data{x,PCA.conIDX}')',[-1 2]);
%     title (['PC: ' num2str(x)]);
% end
% display(['Saved RLD10 to condition sortStruct.idx: ' num2str(PCA.conIDX)]);
% PCA.conIDX=PCA.conIDX+1;