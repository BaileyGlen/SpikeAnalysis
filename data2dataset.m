function [ dataStruct ] = data2dataset( dataStruct )
%DATA2DATASET Takes in the data Struct, and data with dataset
%   Detailed explanation goes here
tempStruct = dataStruct.Output.DPrime.Line;

datasetInc = createTempDataset(tempStruct,'Inc');
datasetDec = createTempDataset(tempStruct,'Dec');
dataStruct.Output.LineDataset =sortrows( vertcat(datasetInc,datasetDec),'Schedule');
%dataStruct.Output.LineDataset=join(datasetInc,datasetDec,{
    function tempDataset = createTempDataset(data,Dir)
        dataFieldNames = fieldnames(dataStruct.Output.DPrime.Line);
        matches = strfind(dataFieldNames,Dir);
        matchingFieldIdx =  ~cellfun(@isempty,matches);
        data = rmfield(data,dataFieldNames(matchingFieldIdx));
        dataFieldNames = dataFieldNames(matchingFieldIdx);
        colNames = cellfun(@(x) x(1:end-8),dataFieldNames,'UniformOutput', false);
        tempDataset = struct2dataset(data);
        tempDataset.Properties.VarNames = colNames;
        if strcmp(Dir,'Inc')
            tempDataset.Direction(1:length(tempDataset),1) = {'increase'};
        elseif strcmp(Dir,'Dec')
            tempDataset.Direction(1:length(tempDataset),1) = {'decrease'};
        else error('The Direction was invalid');
        end
        % setup the session variable
        % posible variables RIRR+day, Deval+Schedule
        % 
        if any(strcmp({'Deval','ContDeg'},data.SessionType))
            
        end
        tempDataset = tempDataset(:,[3 5 4]);
    end
% VarNames = {'CellID','Schedule','Category','zScoreMV'};
% tempDataset = dataset ([],[],[],{},'VarNames',VarNames);
% newDataStruct=tempDataset;

end

