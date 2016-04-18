function [ dataStruct ] = DriverFunction_RRRI_rfdel_z13s_press(  )
%FilesAndSettings This is used to initialize all the various settings, in order to
%have a single place to control them from
%   This is a basic function to set up all the settings for feeding into
%   the other functions


%% Settings for creating the perievent mtx

dataStruct.DataSet='rfdel_z13s_press';
dataStruct.SessionType = 'RRRI';
eventTypeMkPeri='R'; 
eventTypeDPrime='P2';
preEvt=3;            % time prior to Event in sec
postEvt=13;           % time post Event in sec
rasterBin=.1;      % Size of bines for raster
scheduleList={'LL','RL'};
timepointList = {'04', '10'};
if ispc 
    userdir= getenv('USERPROFILE');
    datadir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\_mat\'];
    outputdir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\' dataStruct.DataSet '\'];
else
    userdir= getenv('HOME');
    datadir = '/home/bailey/Documents/MATLAB/RRRI/_mat';
    outputdir = ['/home/bailey/Documents/MATLAB/RRRI/' dataStruct.DataSet];
    %datadir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/_mat/'];
    %outputdir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/z3s_press/'];
end
%datadir='/home/bailey/Documents/MATLAB/RRRI';

%% Create the data Structure
dataStruct.xA=[-1*preEvt:rasterBin:postEvt];

%% Main Loop
dataStruct.fileNameList = {};
for timepointIDX = 1:2
    dataStruct.fileNameList{timepointIDX} = genFileNameList;
    %timepoint = timepointList{timepointIDX}; 
    for scheduleIDX=1:2
        curField=[scheduleList{scheduleIDX} timepointList{timepointIDX}];
        bEvt=[scheduleList{scheduleIDX} '_' eventTypeMkPeri];
        dataStruct.(curField)=getPeriEvent ...
            (bEvt, dataStruct.fileNameList{timepointIDX}, preEvt, postEvt, rasterBin);
        %dataStruct.(curField)=getPCA(dataStruct.(curField));
    end
end
dataStruct=DPrime_2Conditions(dataStruct ,eventTypeDPrime,{'LL04', 'RL04', 'LL10', 'RL10'});
cd (outputdir);
save( dataStruct.DataSet, 'dataStruct');
export(dataStruct.output.dataSet, 'File', [dataStruct.DataSet '.csv'], 'Delimiter', ',');
%% Helper Functions
    function fileNameList = genFileNameList()
        cd (datadir);
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