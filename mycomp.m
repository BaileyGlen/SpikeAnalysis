figure;
subplot(2,2,1);
[sorted,idx] = sort(bsxfun(@rdivide,mean(PCA.LL04.MlMtx(:,10:30),2)-mean(PCA.LL04.MlMtx(:,31:50),2),mean(PCA.LL04.BlMtx,2)));
imagesc(PCA.xA,1:size(sorted,1),flipud(zscore(PCA.LL04.MlMtx(idx,:)')'),[-2 2]);
PCA.LL04.mySort.idx=idx;
subplot(2,2,2);
[sorted,idx] = sort(bsxfun(@rdivide,mean(PCA.RL04.MlMtx(:,10:30),2)-mean(PCA.RL04.MlMtx(:,31:50),2),mean(PCA.RL04.BlMtx,2)));
imagesc(PCA.xA,1:size(sorted,1),flipud(zscore(PCA.RL04.MlMtx(idx,:)')'),[-2 2]);
PCA.RL04.mySort.idx=idx;
subplot(2,2,3);
[sorted,idx] = sort(bsxfun(@rdivide,mean(PCA.LL10.MlMtx(:,10:30),2)-mean(PCA.LL10.MlMtx(:,31:50),2),mean(PCA.LL10.BlMtx,2)));
imagesc(PCA.xA,1:size(sorted,1),flipud(zscore(PCA.LL10.MlMtx(idx,:)')'),[-2 2]);
PCA.LL10.mySort.idx=idx;
subplot(2,2,4);
[sorted,idx] = sort(bsxfun(@rdivide,mean(PCA.RL10.MlMtx(:,10:30),2)-mean(PCA.RL10.MlMtx(:,31:50),2),mean(PCA.RL10.BlMtx,2)));
imagesc(PCA.xA,1:size(sorted,1),flipud(zscore(PCA.RL10.MlMtx(idx,:)')'),[-2 2]);
PCA.RL10.mySort.idx=idx;
%% my comp in reverse
figure;
subplot(2,2,1);
[sorted,idx] = sort(bsxfun(@rdivide,mean(PCA.LL04.MlMtx(:,10:30),2)-mean(PCA.LL04.MlMtx(:,31:50),2),mean(PCA.LL04.BlMtx,2)));
imagesc(PCA.xA,1:size(sorted,1),flipud(zscore(PCA.RL04.MlMtx(idx,:)')'),[-2 2]);
PCA.LL04.mySort.idx=idx;
subplot(2,2,2);
[sorted,idx] = sort(bsxfun(@rdivide,mean(PCA.RL04.MlMtx(:,10:30),2)-mean(PCA.RL04.MlMtx(:,31:50),2),mean(PCA.RL04.BlMtx,2)));
imagesc(PCA.xA,1:size(sorted,1),flipud(zscore(PCA.LL04.MlMtx(idx,:)')'),[-2 2]);
PCA.RL04.mySort.idx=idx;
subplot(2,2,3);
[sorted,idx] = sort(bsxfun(@rdivide,mean(PCA.LL10.MlMtx(:,10:30),2)-mean(PCA.LL10.MlMtx(:,31:50),2),mean(PCA.LL10.BlMtx,2)));
imagesc(PCA.xA,1:size(sorted,1),flipud(zscore(PCA.RL10.MlMtx(idx,:)')'),[-2 2]);
PCA.LL10.mySort.idx=idx;
subplot(2,2,4);
[sorted,idx] = sort(bsxfun(@rdivide,mean(PCA.RL10.MlMtx(:,10:30),2)-mean(PCA.RL10.MlMtx(:,31:50),2),mean(PCA.RL10.BlMtx,2)));
imagesc(PCA.xA,1:size(sorted,1),flipud(zscore(PCA.LL10.MlMtx(idx,:)')'),[-2 2]);
PCA.RL10.mySort.idx=idx;

%%

PCA.LL04.mySort.Inc = mean(zscore(PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:15),:)')');
PCA.RL04.mySort.Inc = mean(zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:15),:)')');
PCA.LL10.mySort.Inc = mean(zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:15),:)')');
PCA.RL10.mySort.Inc = mean(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:15),:)')');

PCA.LL04.mySort.IncErr = std(zscore(PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:15),:)')')/sqrt(15);
PCA.RL04.mySort.IncErr = std(zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:15),:)')')/sqrt(15);
PCA.LL10.mySort.IncErr = std(zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:15),:)')')/sqrt(15);
PCA.RL10.mySort.IncErr = std(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:15),:)')')/sqrt(15);


PCA.LL04.mySort.Dec = mean(zscore(PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-14:end),:)')');
PCA.RL04.mySort.Dec = mean(zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-14:end),:)')');
PCA.LL10.mySort.Dec = mean(zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-14:end),:)')');
PCA.RL10.mySort.Dec = mean(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-14:end),:)')');

PCA.LL04.mySort.DecErr = std(zscore(PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-14:end),:)')')/sqrt(15);
PCA.RL04.mySort.DecErr = std(zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-14:end),:)')')/sqrt(15);
PCA.LL10.mySort.DecErr = std(zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-14:end),:)')')/sqrt(15);
PCA.RL10.mySort.DecErr = std(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-14:end),:)')')/sqrt(15);

%%
PCA.LL04.mySort.Inc = mean((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),:)')');
PCA.RL04.mySort.Inc = mean(zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),:)')');
PCA.LL10.mySort.Inc = mean(zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),:)')');
PCA.RL10.mySort.Inc = mean(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),:)')');

PCA.LL04.mySort.IncErr = std((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),:)')')/sqrt(10);
PCA.RL04.mySort.IncErr = std((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),:)')')/sqrt(10);
PCA.LL10.mySort.IncErr = std((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),:)')')/sqrt(10);
PCA.RL10.mySort.IncErr = std(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),:)')')/sqrt(10);


PCA.LL04.mySort.Dec = mean((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),:)')');
PCA.RL04.mySort.Dec = mean((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),:)')');
PCA.LL10.mySort.Dec = mean((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),:)')');
PCA.RL10.mySort.Dec = mean((PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),:)')');

PCA.LL04.mySort.DecErr = std((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),:)')')/sqrt(10);
PCA.RL04.mySort.DecErr = std((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),:)')')/sqrt(10);
PCA.LL10.mySort.DecErr = std((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),:)')')/sqrt(10);
PCA.RL10.mySort.DecErr = std((PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),:)')')/sqrt(10);

%%
PCA.LL04.mySort.Inc = mean(bsxfun(@rdivide, PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),:), mean(PCA.LL04.BlMtx(PCA.LL04.mySort.idx(1:10),:),2)));
PCA.RL04.mySort.Inc = mean(zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),:)')');
PCA.LL10.mySort.Inc = mean(zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),:)')');
PCA.RL10.mySort.Inc = mean(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),:)')');

PCA.LL04.mySort.IncErr = std((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),:)')')/sqrt(10);
PCA.RL04.mySort.IncErr = std((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),:)')')/sqrt(10);
PCA.LL10.mySort.IncErr = std((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),:)')')/sqrt(10);
PCA.RL10.mySort.IncErr = std(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),:)')')/sqrt(10);


PCA.LL04.mySort.Dec = mean((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),:)')');
PCA.RL04.mySort.Dec = mean((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),:)')');
PCA.LL10.mySort.Dec = mean((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),:)')');
PCA.RL10.mySort.Dec = mean((PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),:)')');

PCA.LL04.mySort.DecErr = std((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),:)')')/sqrt(10);
PCA.RL04.mySort.DecErr = std((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),:)')')/sqrt(10);
PCA.LL10.mySort.DecErr = std((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),:)')')/sqrt(10);
PCA.RL10.mySort.DecErr = std((PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),:)')')/sqrt(10);



%% means and std for pre 1.7 and post 1.7 pretty much useless, see below
preRange = [13:29];
postRange=[30:46];
PCA.LL04.mySort.IncPre = mean(mean(zscore(PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),preRange)')'));
PCA.RL04.mySort.IncPre = mean(mean(zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),preRange)')'));
PCA.LL10.mySort.IncPre = mean(mean(zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),preRange)')'));
PCA.RL10.mySort.IncPre = mean(mean(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),preRange)')'));

PCA.LL04.mySort.IncPost = mean(mean((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),postRange)')'));
PCA.RL04.mySort.IncPost = mean(mean(zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),postRange)')'));
PCA.LL10.mySort.IncPost = mean(mean(zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),postRange)')'));
PCA.RL10.mySort.IncPost = mean(mean(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),postRange)')'));

PCA.LL04.mySort.IncPreErr = mean(std((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),preRange)')')/sqrt(10));
PCA.RL04.mySort.IncPreErr = mean(std((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),preRange)')')/sqrt(10));
PCA.LL10.mySort.IncPreErr = mean(std((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),preRange)')')/sqrt(10));
PCA.RL10.mySort.IncPreErr = mean(std(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),preRange)')')/sqrt(10));

PCA.LL04.mySort.IncPostErr = mean(std((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),postRange)')')/sqrt(10));
PCA.RL04.mySort.IncPostErr = mean(std((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),postRange)')')/sqrt(10));
PCA.LL10.mySort.IncPostErr = mean(std((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),postRange)')')/sqrt(10));
PCA.RL10.mySort.IncPostErr = mean(std(zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),postRange)')')/sqrt(10));

PCA.LL04.mySort.DecPre = mean(mean((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),preRange)')'));
PCA.RL04.mySort.DecPre = mean(mean((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),preRange)')'));
PCA.LL10.mySort.DecPre = mean(mean((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),preRange)')'));
PCA.RL10.mySort.DecPre = mean(mean((PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),preRange)')'));

PCA.LL04.mySort.DecPost = mean(mean((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),postRange)')'));
PCA.RL04.mySort.DecPost = mean(mean((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),postRange)')'));
PCA.LL10.mySort.DecPost = mean(mean((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),postRange)')'));
PCA.RL10.mySort.DecPost = mean(mean((PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),postRange)')'));

PCA.LL04.mySort.DecPreErr = mean(std((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),preRange)')')/sqrt(10));
PCA.RL04.mySort.DecPreErr = mean(std((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),preRange)')')/sqrt(10));
PCA.LL10.mySort.DecPreErr = mean(std((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),preRange)')')/sqrt(10));
PCA.RL10.mySort.DecPreErr = mean(std((PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),preRange)')')/sqrt(10));

PCA.LL04.mySort.DecPostErr = mean(std((PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),postRange)')')/sqrt(10));
PCA.RL04.mySort.DecPostErr = mean(std((PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),postRange)')')/sqrt(10));
PCA.LL10.mySort.DecPostErr = mean(std((PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),postRange)')')/sqrt(10));
PCA.RL10.mySort.DecPostErr = mean(std((PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),postRange)')')/sqrt(10));

var1={'LL04', 'RL04', 'LL10', 'RL10'};
var2={'IncPre', 'IncPost', 'IncPreErr', 'IncPostErr'}
%% second attempt at getting pre and post press means and std in zscore 
IncRaw(1,:,:) = zscore(PCA.LL04.MlMtx(PCA.LL04.mySort.idx(1:10),:)')';
IncRaw(2,:,:) = zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(1:10),:)')';
IncRaw(3,:,:) = zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(1:10),:)')';
IncRaw(4,:,:) = zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(1:10),:)')';

DecRaw(1,:,:) = zscore(PCA.LL04.MlMtx(PCA.LL04.mySort.idx(end-9:end),:)')';
DecRaw(2,:,:) = zscore(PCA.RL04.MlMtx(PCA.RL04.mySort.idx(end-9:end),:)')';
DecRaw(3,:,:) = zscore(PCA.LL10.MlMtx(PCA.LL10.mySort.idx(end-9:end),:)')';
DecRaw(4,:,:) = zscore(PCA.RL10.MlMtx(PCA.RL10.mySort.idx(end-9:end),:)')';
%%
preRange = [13:29];
postRange=[30:46];
for x = 1:4    
    IncRawPre(x,:)=mean(IncRaw(x,:,preRange),3);
    %IncRawSum(x+4,:)=std(IncRaw(x,:,:))/sqrt(10);
    DecRawPre(x,:)=mean(DecRaw(x,:,preRange),3);
    %DecRawSum(x+4,:)=std(DecRaw(x,:,:))/sqrt(10);
    
    IncRawPost(x,:)=mean(IncRaw(x,:,postRange),3);
    %IncRawSum(x+4,:)=std(IncRaw(x,:,:))/sqrt(10);
    DecRawPost(x,:)=mean(DecRaw(x,:,postRange),3);
end
%%
preRange = [13:29];
postRange=[30:46];
IncPre(1:4,:) = mean(IncRawPre(),2);
IncPre(5:8,:) = std(IncRawPre')'/sqrt(10);
IncPost(1:4,:) = mean(IncRawPost(),2);
IncPost(5:8,:) = std(IncRawPost')'/sqrt(10);
DecPre(1:4,:) = mean(DecRawPre(),2);
DecPre(5:8,:) = std(DecRawPre')'/sqrt(10);
DecPost(1:4,:) = mean(DecRawPost(),2);
DecPost(5:8,:) = std(DecRawPost')'/sqrt(10);
