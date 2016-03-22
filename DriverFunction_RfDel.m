function [ dataStruct ] = DriverFunction_RfDel(  )
%FilesAndSettings This is used to initialize all the various settings, in order to
%have a single place to control them from
%   This is a basic function to set up all the settings for feeding into
%   the other functions

%% Settings for creating the perievent mtx
%Generating _04.mat
dataStruct.DataSet='RfDel'; 
eventType='R';
preEvt=3;            % time prior to Event in sec
postEvt=3;           % time post Event in sec
animalList{1}={'m13', 'm14', 'm16','m21', 'm22', 'm25'}; %day 4
animalList{2}={'m13', 'm14', 'm21', 'm22', 'm25'};        %day 10 
fileEnd = '_03.mat';
rasterBin=.1;      % Size of bines for raster

% nmPC=4;
tempRange = [1:((preEvt+postEvt)/rasterBin)+1];
dayList={'04','10'};
scheduleList={'LL','RL'};
mnSpk = 0;
mxSpk=1000;
PLOT=0;
% DayVar={4,4,10,10};
%EventVar={'LL_I','RL_I','LL_I','RL_I'};
%StructVar={'LL04','RL04','LL10','RL10'};
%% removing the above code to instead make this based on...
% animal list, day list, scedule list. build struct var from these



%DirList={'m13D04_03.mat';'m14D04_03.mat';'m16D04_03.mat';'m21D04_03.mat';'m22D04_03.mat';'m25D04_03.mat';...
%            'm13D10_03.mat';'m14D04_03.mat';'m21D10_03.mat';'m22D10_03.mat';'m25D10_03.mat'};



%% Create the data Structure - Reinforced
dataStruct.xA=[-1*preEvt:rasterBin:postEvt];
dataStruct.fileNameList = genFileNameList;
for dayIdx=1:2
    for scheduleIdx=1:2
        curField=[scheduleList{scheduleIdx} dayList{dayIdx}];
        bEvt=[scheduleList{scheduleIdx} '_' eventType];
        dataStruct.(curField)=getPeriEvent;
        dataStruct.(curField)=getPCA(dataStruct.(curField));
    end
end
dataStruct=DPrime_Con(dataStruct ,'P');
% for x=1:4
%     Day=DayVar{x};
%     bEvt = EventVar{x};
%     dataStruct.(StructVar{x})=getPeriEvent;
%     dataStruct.(StructVar{x})=getPCA(dataStruct.(StructVar{x}));
% end
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
    function fileNameList = genFileNameList()
        fileNameList = {};
        for dayIdx=1:2
            for animalIdx=1:length(animalList{dayIdx});
                fileNameList{dayIdx}{animalIdx}=[animalList{dayIdx}{animalIdx} 'D' dayList{dayIdx} fileEnd];
            end
        end
    end
    function eventStruct = getPeriEvent()
        [Mtx, MlMtx, BlMtx, pEvt, pEvt_base]=mkPyeperiEvt(bEvt,'r',dataStruct.fileNameList{dayIdx},'zNo',PLOT, preEvt, postEvt, rasterBin);
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

