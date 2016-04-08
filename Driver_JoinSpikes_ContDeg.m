load ('contDegFileList.mat');
numSessions = length(fileListMAT);
for sessionIDX = 1:numSessions
    fileJoiner(fileListTXT{sessionIDX},fileListMAT{sessionIDX});
end
