function [ dataStruct ] = DriverFunction_RRRI_ISO_PressWindow(  )
%FilesAndSettings This is used to initialize all the various settings, in order to
%have a single place to control them from
%   This is a basic function to set up all the settings for feeding into
%   the other functions


%% Settings for creating the perievent mtx

dataStruct.DataSet='RRRI_ISO_PressWindow';
dataStruct.SessionType = 'RRRI';
eventTypeMkPeri='I';
eventTypeDPrime='P2';
preEvt=3;            % time prior to Event in sec
postEvt=3;           % time post Event in sec
rasterBin=.1;      % Size of bines for raster
scheduleList={'LL','RL'};
timepointList = {'04', '10'};
folderLocation='/home/bailey/Documents/MATLAB/RRRI';

%% Create the data Structure
dataStruct.xA=[-1*preEvt:rasterBin:postEvt];

%% Main Loop
dataStruct.fileNameList = {};
for timepointIDX = 1:2
    dataStruct.fileNameList{timepointIDX} = genFileNameList;
    dataStruct.timepoint = timepointList{timepointIDX}; 
    for scheduleIDX=1:2
        curField=[scheduleList{scheduleIDX} timepointList{timepointIDX}];
        bEvt=[scheduleList{scheduleIDX} '_' eventTypeMkPeri];
        dataStruct.(curField)=getPeriEvent ...
            (bEvt, dataStruct.fileNameList{timepointIDX}, preEvt, postEvt, rasterBin);
        %dataStruct.(curField)=getPCA(dataStruct.(curField));
    end
end
dataStruct=DPrime_2Conditions(dataStruct ,eventTypeDPrime,{'LL04', 'RL04', 'LL10', 'RL10'});


%% Helper Functions
    function fileNameList = genFileNameList()
        cd (folderLocation);
        fileNameStruct = dir(['*D' timepointList{timepointIDX} '*.mat']);
        for fileIDX = 1:length(fileNameStruct)
            fileNameList{1,fileIDX} = fileNameStruct(fileIDX).name;
        end
%         fileNameStruct = dir(['*D' timepointList{timepointIDX}]);
%         for fileIDX = 1:length(fileNameStruct)
%             fileNameList{2,fileIDX} = fileNameStruct(fileIDX).name;
%         end
    end
    function eventStruct = getPeriEvent(bEvt, fileNameList, preEvt, postEvt, rasterBin)
        [eventStruct]=mkPeriEvt(bEvt,fileNameList, preEvt, postEvt, rasterBin);
    end
    function eventStruct = getPCA (eventStruct)
        k=find(sum(eventStruct.MlMtx') & sum(eventStruct.MlMtx')); % get rid of low FR ST
        [coeff,score,latent,tsquared,explained,mu]=pca(bsxfun(@rdivide, eventStruct.MlMtx(k,dataStruct.xA),nanmean(eventStruct.BlMtx(k,:),2)));
        %[coeff,score,latent,tsquared,explained,mu]=pca(MlMtx(k,:));
        eventStruct.PCA.coeff=coeff;
        eventStruct.PCA.score=score;
        eventStruct.PCA.latent=latent;
        eventStruct.PCA.tsquared=tsquared;
        eventStruct.PCA.explained=explained;
        eventStruct.PCA.mu=mu;
    end
end