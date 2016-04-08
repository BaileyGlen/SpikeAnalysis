%function [ dataStruct ] = DriverFunction_Deval_ISO(  )
%FilesAndSettings This is used to initialize all the various settings, in order to
%have a single place to control them from
%   This is a basic function to set up all the settings for feeding into
%   the other functions


Load the summary datafile


%% General Settings 

dataStruct.SessionType = 'Deval';
data.leverList =  {'Left','Right'};
data.scheduleList = {'RI','RR'};
folderLocation='/home/bailey/Documents/MATLAB/RRRI';

%% Isolated Presses

dataStruct.DataSet='Deval_ISO';
eventTypeMkPeri='I';
eventTypeDPrime='P';
preEvt=3;            % time prior to Event in sec
postEvt=3;           % time post Event in sec
rasterBin=.1;      % Size of bines for raster


%% Create the data Structure
dataStruct.xA=[-1*preEvt:rasterBin:postEvt];
dataStruct.fileNameList = genFileNameList;
%% Main Loop
for 1:length(fileNameList;
for 
    curField=[scheduleList{scheduleIDX}];
    bEvt=[curField '_' eventTypeMkPeri];
    dataStruct.(curField)=getPeriEvent(bEvt,dataStruct.fileNameList(scheduleIDX,:),preEvt,postEvt,rasterBin);
    %dataStruct.(curField)=getPCA(dataStruct.(curField));
end
dataStruct=DPrime_2Conditions(dataStruct ,eventTypeDPrime);


%% Helper Functions
    function fileNameList = genFileNameList()
        fileNameList{1,:} = getDirectoryFilenames (folderLocation,[],'*RRRI*');
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