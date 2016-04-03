load ('devalFileList.mat');
numSessions = length(fileListMat);
for sessionIDX = 1:numSessions
    fileJoiner(fileListTXT{sessionIDX},fileListMat{sessionIDX});
end