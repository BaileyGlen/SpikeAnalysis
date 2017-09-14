function [ dataStruct ] = DriverFunction_RRRI_hz(  )
%FilesAndSettings This is used to initialize all the various settings, in order to
%have a single place to control them from
%   This is a basic function to set up all the settings for feeding into
%   the other functions


%% Settings for creating the perievent mtx

dataStruct.DataSet='RRRI';
timepointList = {'04', '10'};
if ispc 
    userdir= getenv('USERPROFILE');
    datadir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\_mat\'];
    %outputdir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\' dataStruct.DataSet '\'];
else
    userdir= getenv('HOME');
    datadir = '/home/bailey/Documents/MATLAB/RRRI/_mat';
    %outputdir = ['/home/bailey/Documents/MATLAB/RRRI/' dataStruct.DataSet];
    %datadir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/_mat/'];
    %outputdir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/z3s_press/'];
end
%datadir='/home/bailey/Documents/MATLAB/RRRI';

%% Create the data Structure
%% Main Loop
dataStruct.fileNameList = {};
datasetHZ = dataset()
for timepointIDX = 1:2
    dataStruct.fileNameList{timepointIDX} = genFileNameList;
    %timepoint = timepointList{timepointIDX}; 
    datasetHZ = getHZ(dataStruct.fileNameList{timepointIDX}, datasetHZ, timepointList{timepointIDX});
end

cd ('C:\Users\Jacqui\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\HZ');
save( 'datasetHZ.mat', 'datasetHZ');
%export(dataStruct.output.dataSet, 'File', [dataStruct.DataSet '.csv'], 'Delimiter', ',');

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

 
end
