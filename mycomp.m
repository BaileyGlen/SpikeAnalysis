figure;
subplot(2,2,1);
[sorted,idx] = sort(bsxfun(@rdivide,mean(data.LL04.MlMtx(:,10:30),2)-mean(data.LL04.MlMtx(:,31:50),2),mean(data.LL04.BlMtx,2)));
imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.LL04.MlMtx(idx,:)')'),[-2 2]);
data.LL04.mySort.idx=idx;
subplot(2,2,2);
[sorted,idx] = sort(bsxfun(@rdivide,mean(data.RL04.MlMtx(:,10:30),2)-mean(data.RL04.MlMtx(:,31:50),2),mean(data.RL04.BlMtx,2)));
imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.RL04.MlMtx(idx,:)')'),[-2 2]);
data.RL04.mySort.idx=idx;
subplot(2,2,3);
[sorted,idx] = sort(bsxfun(@rdivide,mean(data.LL10.MlMtx(:,10:30),2)-mean(data.LL10.MlMtx(:,31:50),2),mean(data.LL10.BlMtx,2)));
imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.LL10.MlMtx(idx,:)')'),[-2 2]);
data.LL10.mySort.idx=idx;
subplot(2,2,4);
[sorted,idx] = sort(bsxfun(@rdivide,mean(data.RL10.MlMtx(:,10:30),2)-mean(data.RL10.MlMtx(:,31:50),2),mean(data.RL10.BlMtx,2)));
imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.RL10.MlMtx(idx,:)')'),[-2 2]);
data.RL10.mySort.idx=idx;
%% my comp in reverse
figure;
subplot(2,2,1);
[sorted,idx] = sort(bsxfun(@rdivide,mean(data.LL04.MlMtx(:,10:30),2)-mean(data.LL04.MlMtx(:,31:50),2),mean(data.LL04.BlMtx,2)));
imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.RL04.MlMtx(idx,:)')'),[-2 2]);
data.LL04.mySort.idx=idx;
subplot(2,2,2);
[sorted,idx] = sort(bsxfun(@rdivide,mean(data.RL04.MlMtx(:,10:30),2)-mean(data.RL04.MlMtx(:,31:50),2),mean(data.RL04.BlMtx,2)));
imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.LL04.MlMtx(idx,:)')'),[-2 2]);
data.RL04.mySort.idx=idx;
subplot(2,2,3);
[sorted,idx] = sort(bsxfun(@rdivide,mean(data.LL10.MlMtx(:,10:30),2)-mean(data.LL10.MlMtx(:,31:50),2),mean(data.LL10.BlMtx,2)));
imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.RL10.MlMtx(idx,:)')'),[-2 2]);
data.LL10.mySort.idx=idx;
subplot(2,2,4);
[sorted,idx] = sort(bsxfun(@rdivide,mean(data.RL10.MlMtx(:,10:30),2)-mean(data.RL10.MlMtx(:,31:50),2),mean(data.RL10.BlMtx,2)));
imagesc(data.xA,1:size(sorted,1),flipud(zscore(data.LL10.MlMtx(idx,:)')'),[-2 2]);
data.RL10.mySort.idx=idx;

%%

data.LL04.mySort.Inc = mean(zscore(data.LL04.MlMtx(data.LL04.mySort.idx(1:15),:)')');
data.RL04.mySort.Inc = mean(zscore(data.RL04.MlMtx(data.RL04.mySort.idx(1:15),:)')');
data.LL10.mySort.Inc = mean(zscore(data.LL10.MlMtx(data.LL10.mySort.idx(1:15),:)')');
data.RL10.mySort.Inc = mean(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:15),:)')');

data.LL04.mySort.IncErr = std(zscore(data.LL04.MlMtx(data.LL04.mySort.idx(1:15),:)')')/sqrt(15);
data.RL04.mySort.IncErr = std(zscore(data.RL04.MlMtx(data.RL04.mySort.idx(1:15),:)')')/sqrt(15);
data.LL10.mySort.IncErr = std(zscore(data.LL10.MlMtx(data.LL10.mySort.idx(1:15),:)')')/sqrt(15);
data.RL10.mySort.IncErr = std(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:15),:)')')/sqrt(15);


data.LL04.mySort.Dec = mean(zscore(data.LL04.MlMtx(data.LL04.mySort.idx(end-14:end),:)')');
data.RL04.mySort.Dec = mean(zscore(data.RL04.MlMtx(data.RL04.mySort.idx(end-14:end),:)')');
data.LL10.mySort.Dec = mean(zscore(data.LL10.MlMtx(data.LL10.mySort.idx(end-14:end),:)')');
data.RL10.mySort.Dec = mean(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(end-14:end),:)')');

data.LL04.mySort.DecErr = std(zscore(data.LL04.MlMtx(data.LL04.mySort.idx(end-14:end),:)')')/sqrt(15);
data.RL04.mySort.DecErr = std(zscore(data.RL04.MlMtx(data.RL04.mySort.idx(end-14:end),:)')')/sqrt(15);
data.LL10.mySort.DecErr = std(zscore(data.LL10.MlMtx(data.LL10.mySort.idx(end-14:end),:)')')/sqrt(15);
data.RL10.mySort.DecErr = std(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(end-14:end),:)')')/sqrt(15);

%%
data.LL04.mySort.Inc = mean((data.LL04.MlMtx(data.LL04.mySort.idx(1:10),:)')');
data.RL04.mySort.Inc = mean(zscore(data.RL04.MlMtx(data.RL04.mySort.idx(1:10),:)')');
data.LL10.mySort.Inc = mean(zscore(data.LL10.MlMtx(data.LL10.mySort.idx(1:10),:)')');
data.RL10.mySort.Inc = mean(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),:)')');

data.LL04.mySort.IncErr = std((data.LL04.MlMtx(data.LL04.mySort.idx(1:10),:)')')/sqrt(10);
data.RL04.mySort.IncErr = std((data.RL04.MlMtx(data.RL04.mySort.idx(1:10),:)')')/sqrt(10);
data.LL10.mySort.IncErr = std((data.LL10.MlMtx(data.LL10.mySort.idx(1:10),:)')')/sqrt(10);
data.RL10.mySort.IncErr = std(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),:)')')/sqrt(10);


data.LL04.mySort.Dec = mean((data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),:)')');
data.RL04.mySort.Dec = mean((data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),:)')');
data.LL10.mySort.Dec = mean((data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),:)')');
data.RL10.mySort.Dec = mean((data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),:)')');

data.LL04.mySort.DecErr = std((data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),:)')')/sqrt(10);
data.RL04.mySort.DecErr = std((data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),:)')')/sqrt(10);
data.LL10.mySort.DecErr = std((data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),:)')')/sqrt(10);
data.RL10.mySort.DecErr = std((data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),:)')')/sqrt(10);

%%
data.LL04.mySort.Inc = mean(bsxfun(@rdivide, data.LL04.MlMtx(data.LL04.mySort.idx(1:10),:), mean(data.LL04.BlMtx(data.LL04.mySort.idx(1:10),:),2)));
data.RL04.mySort.Inc = mean(zscore(data.RL04.MlMtx(data.RL04.mySort.idx(1:10),:)')');
data.LL10.mySort.Inc = mean(zscore(data.LL10.MlMtx(data.LL10.mySort.idx(1:10),:)')');
data.RL10.mySort.Inc = mean(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),:)')');

data.LL04.mySort.IncErr = std((data.LL04.MlMtx(data.LL04.mySort.idx(1:10),:)')')/sqrt(10);
data.RL04.mySort.IncErr = std((data.RL04.MlMtx(data.RL04.mySort.idx(1:10),:)')')/sqrt(10);
data.LL10.mySort.IncErr = std((data.LL10.MlMtx(data.LL10.mySort.idx(1:10),:)')')/sqrt(10);
data.RL10.mySort.IncErr = std(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),:)')')/sqrt(10);


data.LL04.mySort.Dec = mean((data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),:)')');
data.RL04.mySort.Dec = mean((data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),:)')');
data.LL10.mySort.Dec = mean((data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),:)')');
data.RL10.mySort.Dec = mean((data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),:)')');

data.LL04.mySort.DecErr = std((data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),:)')')/sqrt(10);
data.RL04.mySort.DecErr = std((data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),:)')')/sqrt(10);
data.LL10.mySort.DecErr = std((data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),:)')')/sqrt(10);
data.RL10.mySort.DecErr = std((data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),:)')')/sqrt(10);



%% means and std for pre 1.7 and post 1.7 pretty much useless, see below
preRange = [13:29];
postRange=[30:46];
data.LL04.mySort.IncPre = mean(mean(zscore(data.LL04.MlMtx(data.LL04.mySort.idx(1:10),preRange)')'));
data.RL04.mySort.IncPre = mean(mean(zscore(data.RL04.MlMtx(data.RL04.mySort.idx(1:10),preRange)')'));
data.LL10.mySort.IncPre = mean(mean(zscore(data.LL10.MlMtx(data.LL10.mySort.idx(1:10),preRange)')'));
data.RL10.mySort.IncPre = mean(mean(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),preRange)')'));

data.LL04.mySort.IncPost = mean(mean((data.LL04.MlMtx(data.LL04.mySort.idx(1:10),postRange)')'));
data.RL04.mySort.IncPost = mean(mean(zscore(data.RL04.MlMtx(data.RL04.mySort.idx(1:10),postRange)')'));
data.LL10.mySort.IncPost = mean(mean(zscore(data.LL10.MlMtx(data.LL10.mySort.idx(1:10),postRange)')'));
data.RL10.mySort.IncPost = mean(mean(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),postRange)')'));

data.LL04.mySort.IncPreErr = mean(std((data.LL04.MlMtx(data.LL04.mySort.idx(1:10),preRange)')')/sqrt(10));
data.RL04.mySort.IncPreErr = mean(std((data.RL04.MlMtx(data.RL04.mySort.idx(1:10),preRange)')')/sqrt(10));
data.LL10.mySort.IncPreErr = mean(std((data.LL10.MlMtx(data.LL10.mySort.idx(1:10),preRange)')')/sqrt(10));
data.RL10.mySort.IncPreErr = mean(std(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),preRange)')')/sqrt(10));

data.LL04.mySort.IncPostErr = mean(std((data.LL04.MlMtx(data.LL04.mySort.idx(1:10),postRange)')')/sqrt(10));
data.RL04.mySort.IncPostErr = mean(std((data.RL04.MlMtx(data.RL04.mySort.idx(1:10),postRange)')')/sqrt(10));
data.LL10.mySort.IncPostErr = mean(std((data.LL10.MlMtx(data.LL10.mySort.idx(1:10),postRange)')')/sqrt(10));
data.RL10.mySort.IncPostErr = mean(std(zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),postRange)')')/sqrt(10));

data.LL04.mySort.DecPre = mean(mean((data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),preRange)')'));
data.RL04.mySort.DecPre = mean(mean((data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),preRange)')'));
data.LL10.mySort.DecPre = mean(mean((data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),preRange)')'));
data.RL10.mySort.DecPre = mean(mean((data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),preRange)')'));

data.LL04.mySort.DecPost = mean(mean((data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),postRange)')'));
data.RL04.mySort.DecPost = mean(mean((data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),postRange)')'));
data.LL10.mySort.DecPost = mean(mean((data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),postRange)')'));
data.RL10.mySort.DecPost = mean(mean((data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),postRange)')'));

data.LL04.mySort.DecPreErr = mean(std((data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),preRange)')')/sqrt(10));
data.RL04.mySort.DecPreErr = mean(std((data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),preRange)')')/sqrt(10));
data.LL10.mySort.DecPreErr = mean(std((data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),preRange)')')/sqrt(10));
data.RL10.mySort.DecPreErr = mean(std((data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),preRange)')')/sqrt(10));

data.LL04.mySort.DecPostErr = mean(std((data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),postRange)')')/sqrt(10));
data.RL04.mySort.DecPostErr = mean(std((data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),postRange)')')/sqrt(10));
data.LL10.mySort.DecPostErr = mean(std((data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),postRange)')')/sqrt(10));
data.RL10.mySort.DecPostErr = mean(std((data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),postRange)')')/sqrt(10));

var1={'LL04', 'RL04', 'LL10', 'RL10'};
var2={'IncPre', 'IncPost', 'IncPreErr', 'IncPostErr'}
%% second attempt at getting pre and post press means and std in zscore 
IncRaw(1,:,:) = zscore(data.LL04.MlMtx(data.LL04.mySort.idx(1:10),:)')';
IncRaw(2,:,:) = zscore(data.RL04.MlMtx(data.RL04.mySort.idx(1:10),:)')';
IncRaw(3,:,:) = zscore(data.LL10.MlMtx(data.LL10.mySort.idx(1:10),:)')';
IncRaw(4,:,:) = zscore(data.RL10.MlMtx(data.RL10.mySort.idx(1:10),:)')';

DecRaw(1,:,:) = zscore(data.LL04.MlMtx(data.LL04.mySort.idx(end-9:end),:)')';
DecRaw(2,:,:) = zscore(data.RL04.MlMtx(data.RL04.mySort.idx(end-9:end),:)')';
DecRaw(3,:,:) = zscore(data.LL10.MlMtx(data.LL10.mySort.idx(end-9:end),:)')';
DecRaw(4,:,:) = zscore(data.RL10.MlMtx(data.RL10.mySort.idx(end-9:end),:)')';
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
