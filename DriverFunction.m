function [ dataStruct ] = DriverFunction(  )
%FilesAndSettings This is used to initialize all the various settings, in order to
%have a single place to control them from
%   This is a basic function to set up all the settings for feeding into
%   the other functions

%% Settings for creating the perievent mtx
mnSpk = 0;
mxSpk=1000;
rasterBin=.1;      % Size of bines for raster
preEvt=3;            % time prior to Event in sec
postEvt=10;           % time post Event in sec
PLOT=0;
nmPC=4;
tempRange = [1:((preEvt+postEvt)/rasterBin)+1];
DayVar={4,4,10,10};
EventVar={'LL_R','RL_R','LL_R','RL_R'};
StructVar={'LL04','RL04','LL10','RL10'};
%% Create the data Structure - Reinforced
dataStruct.xA=[-1*preEvt:rasterBin:postEvt];
for x=1:4
    Day=DayVar{x};
    bEvt = EventVar{x};
    dataStruct.(StructVar{x})=getPeriEvent;
    dataStruct.(StructVar{x})=getPCA(dataStruct.(StructVar{x}));
end
% Day=4;
% bEvt1={'LL_R'};
% dataStruct.LL04=getPeriEvent;
% dataStruct.LL04=getPCA;
% Day=4;
% bEvt1={'RL_R'};
% dataStruct.RL04=getPeriEvent;
% dataStruct.RL04=getPCA;
% Day=10;
% bEvt1={'RL_R'};
% dataStruct.LL10=getPeriEvent;
% dataStruct.LL10=getPCA;
% Day=10;
% bEvt1={'RL_R'};
% dataStruct.RL10=getPeriEvent;
% dataStruct.RL10=getPCA;


    function eventStruct = getPeriEvent()
        [Mtx, MlMtx, BlMtx, pEvt, pEvt_base]=mkPeriEvt(bEvt,'r',Day,'zNo',PLOT, preEvt, postEvt, rasterBin);
        eventStruct.Mtx = Mtx;
        eventStruct.MlMtx = MlMtx;
        eventStruct.BlMtx = BlMtx;
        eventStruct.pEvt = pEvt;
        eventStruct.pEvt_base = pEvt_base;
    end
    function eventStruct = getPCA (eventStruct)
        k=find(sum(eventStruct.MlMtx')>mnSpk & sum(eventStruct.MlMtx')<mxSpk); % get rid of low FR ST
        [coeff,score,latent,tsquared,explained,mu]=pca(bsxfun(@rdivide, eventStruct.MlMtx(k,tempRange),nanmean(eventStruct.BlMtx(k,:),2)));
        %[coeff,score,latent,tsquared,explained,mu]=pca(MlMtx(k,:));
        eventStruct.PCA.coeff=coeff;
        eventStruct.PCA.score=score;
        eventStruct.PCA.latent=latent;
        eventStruct.PCA.tsquared=tsquared;
        eventStruct.PCA.explained=explained;
        eventStruct.PCA.mu=mu;
        
    end

end

