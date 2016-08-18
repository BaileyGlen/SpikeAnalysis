function [ dataStruct ] = DriverFunction_ContDeg_RfDel(  )
%FilesAndSettings This is used to initialize all the various settings, in order to
%have a single place to control them from
%   This is a basic function to set up all the settings for feeding into
%   the other functions


%% Settings for creating the perievent mtx

dataStruct.DataSet='ContDeg_RfDel';
dataStruct.SessionType = 'ContDeg';
eventTypeMkPeri='Sucrose';
eventTypeDPrime='P';
preEvt=3;            % time prior to Event in sec
postEvt=10;           % time post Event in sec
rasterBin=.1;      % Size of bines for raster
scheduleList={'LL','RL'};
userdir= getenv('USERPROFILE');
    datadir = [userdir '\Box Sync\mea_data\rrri_01\processing\contdeg\spiketrain\_mat\'];
    outputdir = [userdir '\Box Sync\mea_data\rrri_01\processing\contdeg\spiketrain\' dataStruct.DataSet '\'];

%% Create the data Structure
dataStruct.xA=[-1*preEvt:rasterBin:postEvt];
dataStruct.fileNameList = genFileNameList;
%% Main Loop
for scheduleIDX=1:2
    curField=[scheduleList{scheduleIDX}];
    bEvt=[curField '_' eventTypeMkPeri];
    dataStruct.(curField)=getPeriEvent(bEvt,dataStruct.fileNameList(scheduleIDX,:),preEvt,postEvt,rasterBin);
    %dataStruct.(curField)=getPCA(dataStruct.(curField));
end
dataStruct=DPrime_2Conditions(dataStruct ,eventTypeDPrime, {'LL', 'RL'});
cd (outputdir);
save( dataStruct.DataSet, 'dataStruct');
export(dataStruct.output.dataSet, 'File', [dataStruct.DataSet '.csv'], 'Delimiter', ',');


%% Helper Functions
    function fileNameList = genFileNameList()
        cd (datadir);
        fileNameStruct = dir('*RI*final*');
        for fileIDX = 1:length(fileNameStruct)
            fileNameList{1,fileIDX} = fileNameStruct(fileIDX).name;
        end
        fileNameStruct = dir('*RR*final*');
        for fileIDX = 1:length(fileNameStruct)
            fileNameList{2,fileIDX} = fileNameStruct(fileIDX).name;
        end
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

