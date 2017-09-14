function [ dataStruct ] = DriverFunction_Deval_hz(  )
%FilesAndSettings This is used to initialize all the various settings, in order to
%have a single place to control them from
%   This is a basic function to set up all the settings for feeding into
%   the other functions


%% Settings for creating the perievent mtx

dataStruct.DataSet='Deval';
scheduleList={'LL','RL'};
if ispc 
    userdir= getenv('USERPROFILE');
    datadir = [userdir '\Box Sync\mea_data\rrri_01\processing\deval\spiketrain\_mat\'];
    %outputdir = [userdir '\Box Sync\mea_data\rrri_01\processing\rrri\spiketrain\' dataStruct.DataSet '\'];
else
    userdir= getenv('HOME');
    datadir = '/home/bailey/Documents/MATLAB/deval/_mat';
    %outputdir = ['/home/bailey/Documents/MATLAB/RRRI/' dataStruct.DataSet];
    %datadir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/_mat/'];
    %outputdir = [userdir '/Box Sync/mea_data/rrri_01/processing/rrri/spiketrain/z3s_press/'];
end
%datadir='/home/bailey/Documents/MATLAB/RRRI';

%% Create the data Structure
%% Main Loop
dataStruct.fileNameList = {};
datasetHZ = dataset()
dataStruct.fileNameList = genFileNameList;
for scheduleIDX = 1:2
    %timepoint = timepointList{timepointIDX}; 
    datasetHZ = getHZ_testSession(dataStruct.fileNameList(scheduleIDX,:), datasetHZ, scheduleList{scheduleIDX});
end

cd ('C:\Users\Jacqui\Box Sync\mea_data\rrri_01\processing\deval\spiketrain\HZ');
save( 'datasetHZ.mat', 'datasetHZ');
%export(dataStruct.output.dataSet, 'File', [dataStruct.DataSet '.csv'], 'Delimiter', ',');

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
 
end
