% this is to run the hmmTimeSegmentation code
INITIAL_NUM_STATES = 6;
DIM = 19;
MAX_NUM_FEATURES = 500000;
allData = zeros(MAX_NUM_FEATURES, DIM);
allFeatureFlag = zeros(MAX_NUM_FEATURES, 2);
allFeatureFlag(:, 1) = 1;
allStateMeans = zeros(INITIAL_NUM_STATES, DIM);
allStateVars = zeros(INITIAL_NUM_STATES, DIM, DIM);
numStates = INITIAL_NUM_STATES;
MIN_DUR = 200;

%[allData, framePeriod, dataType, typeCode, textVersion] = readhtk('/home/jitu/btp/data/diarization_data/mfcc/IS1001d.fea');
[allData, framePeriod, dataType, typeCode, textVersion] = readhtk('/home/jitu/btp/data/timit_dataset/2.fea');
%==========================================================================
% Classify feature vectors as voiced and unvoiced
allFeatureCount = size(allData);
%[featureCount, data, allFeatureFlag] = ClassificationVoicedUnvoiced(allData, allFeatureFlag);
featureCount = allFeatureCount;
voicedData = zeros(featureCount, DIM);
% for i = 1:featureCount
%     voicedData(i, :) = data(i, :);
% end
for i = 1:featureCount
    voicedData(i, :) = allData(i, :);
end
% regularization inverse cov. matrices
%function [hmm,pth,comb] = hmmtimesegment(x,N,MD,stopcrit,reg)
% 
 %Default we use:  STOPCRIT.maxiters = 10
%                  STOPCRIT.minllimpr = 1e-5
%                  REG = 1e-5;
% See also: mogp, moginit
stopcrit.maxiters = 10;
stopcrit.minllimpr = 1e-5;

[hmm, pth, comb] = hmmtimesegment(voicedData, INITIAL_NUM_STATES, MIN_DUR, stopcrit, 1e-3);
disp('experiment done...');
seq = pth;
for i = 1:featureCount
    disp(seq(i, 1));
end
%==========================================================================
% writing the final data in RTTM file in the format specified by NIST
% allData , allFeatureFlag, numStates, fileName
file_name = 'IS1007a';
extension = '.rttm';
rttm_file = fopen(strcat(file_name, extension), 'w');
for i = 1:allFeatureCount
    allFeatureFlag(i, 2) = i;
end

allStateSequence = zeros(allFeatureCount, 1);
for k = 1:allFeatureCount(1,1)
    if allFeatureFlag(k, 2) == 0
        allStateSequence(k, 1) = 0;
    else
        index = allFeatureFlag(k, 2);
        state = path(index, 1);
        allStateSequence(k, 1) = state;
    end
end
for s = 1:numStates
    spkr_name = strcat(file_name, '_', num2str(s));
    fprintf(rttm_file, 'SPKR-INFO %s 1 <NA> <NA> <NA> unknown %s <NA>\n', file_name, spkr_name);
    start = 0; dur = 0;
    for k = 1:allFeatureCount(1, 1)
        disp(k);
        if (allStateSequence(k, 1) ~= s && dur > 0) || k == featureCount
            % write into RTTM file
            fprintf(rttm_file, 'SPEAKER %s 1 %f %f <NA> <NA> %s <NA>\n', file_name, start * 0.010 , dur * 0.010, spkr_name);
            dur = 0;
        end
        
        if allStateSequence(k, 1) == s && dur == 0
            start = k;
        end
        
        if allStateSequence(k, 1) == s
            dur = dur + 1;
        end
    end
end

