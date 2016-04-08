
fileNameStruct = dir('*RI*final*');
for fileIDX = 1:length(fileNameStruct)
    load(fileNameStruct(fileIDX).name);
    
    num2Add = length(data.LL_Rf_TS);
    newTm = [data.Lapish.behaveEvtTm_Raw; data.LL_Rf_TS];
    temp = cell(num2Add,1);
    temp(:) = {'RF'};
    newName = [data.Lapish.behaveEvt_Raw; temp];
    
    [data.Lapish.behaveEvtTm_Raw, sortIDX] = sort(newTm);
    data.Lapish.behaveEvt_Raw= newName(sortIDX);
    save(fileNameStruct(fileIDX).name,'data');
    display(['Processing ' fileNameStruct(fileIDX).name]);
end

%%
fileNameStruct = dir('*RR*final*');
for fileIDX = 1:length(fileNameStruct)
    load(fileNameStruct(fileIDX).name);
    
    num2Add = length(data.RR_Rf_TS);
    newTm = [data.Lapish.behaveEvtTm_Raw; data.RR_Rf_TS];
    temp = cell(num2Add,1);
    temp(:) = {'RF'};
    newName = [data.Lapish.behaveEvt_Raw; temp];
    
    [data.Lapish.behaveEvtTm_Raw, sortIDX] = sort(newTm);
    data.Lapish.behaveEvt_Raw= newName(sortIDX);
    save(fileNameStruct(fileIDX).name,'data');
    display(['Processing ' fileNameStruct(fileIDX).name]);
end
