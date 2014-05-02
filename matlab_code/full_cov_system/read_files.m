% read the HTK format feature file
% function [d,fp,dt,tc,t]=readhtk(file)
%READHTK  read an HTK parameter file [D,FP,DT,TC,T]=(FILE)
INITIAL_NUM_STATES = 20;
DIM = 19;
MAX_NUM_FEATURES = 500000;
allData = zeros(MAX_NUM_FEATURES, DIM);
allFeatureFlag = zeros(MAX_NUM_FEATURES, 2);
allStateMeans = zeros(INITIAL_NUM_STATES, DIM);
allStateVars = zeros(INITIAL_NUM_STATES, DIM, DIM);
numStates = INITIAL_NUM_STATES;

[allData, framePeriod, dataType, typeCode, textVersion] = readhtk('/home/jitu/btp/data/diarization_data/mfcc/IS1007a.fea');
%==========================================================================
% Classify feature vectors as voiced and unvoiced
allFeatureCount = size(allData);
[featureCount, data, allFeatureFlag] = ClassificationVoicedUnvoiced(allData, allFeatureFlag);
%==========================================================================
% Initialize each STATE with a Gaussian having full covariance matrix
% [allStateMeans, allStateVars] = InitializeAllStates(data containing all voiced feature vectors, initial
% no of states, array containg means for each state, matrix containing
% covariance matrix for each state) ;
%[allStateMeans, allStateVars] = InitializeAllStates(data, INITIAL_NUM_STATES, allStateMeans, allStateVars, featureCount);
j = 1;
block = featureCount/ INITIAL_NUM_STATES;
block = int32(block);
x = zeros(block, 19);
stateSequence = zeros(featureCount, 1);
for s=1:INITIAL_NUM_STATES
    i = 1;
    for j=(s-1)*block+1:s*block
        x(i, :) = data(j, :);
        stateSequence(j, 1) = s;
        i = i+1;
    end
    disp(size(x));    
    mu = mean(x);
    c = cov(x);
    % we should store obtained mean and variance into all mixture means and
    % variances
    allStateMeans(s, :) = mu(1, :);
    allStateVars(s, :, :) = c(:, :);
end

voicedData = zeros(featureCount, DIM);
for i = 1:featureCount
    voicedData(i, :) = data(i, :);
end
log_prob = zeros( INITIAL_NUM_STATES, featureCount);

% calculate log probability for every feature vector
for s = 1: numStates
    sigma = squeeze(allStateVars(s, :, :));
    prob = mvnpdf(voicedData, allStateMeans(s, :), sigma);
    prob = log(prob);
    log_prob(s, :) = prob;
end
numElemEachState = zeros(numStates, 1);
pi = zeros(numStates, 1);
for i = 1: numStates
    numElemEachState(i, 1) = int32(featureCount/numStates);
end
%

for i = 1: numStates
    pi(i, 1) = log(numElemEachState(i, 1)/featureCount);
end
%==========================================================================
% Clustering and Merging step :- this step will include viterbi realignment
% , Reclustering of all states using Expectation maximization, calculation
% of deltaBIC, Merging of two STATES
MAX_MERGE_ITER = 6;
mergeIter = 0;
iterFlag = 1;
while numStates > 0 && mergeIter < MAX_MERGE_ITER && iterFlag == 1
    % VITERBI REALIGNMENT
    %[stateSequence, numElemEachState] = viterbi_realign(log_prob, pi, numStates, featureCount, numElemEachState);
    %==========================================================================
    % APPLY MIN_DURATION constraint on viterbi realignment
    % MIN_DUR = 250;
    % for s = 1:numStates
    %     if numElemEachState(s, 1) < MIN_DUR
    %         % we need to drop this state
    %
    %     else
    %     end
    % end
    %==========================================================================
    % reclustering of each state
    disp('reculstering....');
    for s = 1:numStates
        %find all the feature vectors which belong to state s
        count = 1;
        clear tempData;
        for i = 1:featureCount
            if stateSequence(i, 1) == s
                %add it to the tempData
                tempData(count, :) = data(i, :);
                count = count + 1;
            end
        end
        
        %find mean and covariance matrix for this tempData
        mu = mean(tempData);
        sigma = cov(tempData);
        % we should store obtained mean and variance into all mixture means and
        % variances
        allStateMeans(s, :) = mu(1, :);
        allStateVars(s, :, :) = sigma(:, :);
    end
    %==========================================================================
    % calculate delta BIC for each pair of states
    disp('calculating delta BIC for each pair of states');
    bic = zeros(numStates, numStates);
    for i = 1: numStates
        for j = i+1:numStates
            % Model both the states using same gaussian
            count = 1;
            clear tempData;
            for k = 1:featureCount
                if stateSequence(k, 1) == i || stateSequence(k, 1) == j
                    %disp('executed');
                    tempData(count, :) = data(k, :);
                    count = count + 1;
                end
            end
            mu = mean(tempData);
            sigma = cov(tempData);
            sigma_i = squeeze(allStateVars(i, :, :));
            sigma_j = squeeze(allStateVars(j, :, :));
            LAMBDA = 9.0;
            P = 0.5 * (DIM + 0.5 * DIM * (DIM+1)) * log(numElemEachState(i, 1) + numElemEachState(j, 1));
            bic(i, j) = (numElemEachState(i, 1) + numElemEachState(j, 1)) * log(det(sigma)) - (numElemEachState(i,1)) * log(det(sigma_i)) - numElemEachState(j, 1) * log(det(sigma_j)) - LAMBDA * P;
            disp('bic:');
            disp(bic(i, j));
        end
    end
    %==========================================================================
    %find the state with minimum delta BIC
    min_bic = 99999999;
    for i = 1:numStates
        for j = i+1:numStates
            if bic(i, j) < min_bic
                min_bic  = bic(i, j);
                min_i = i;
                min_j = j;
            end
        end
    end
    %==========================================================================
    % stopping criterion  Merge Two States
    flag = 0;
    disp(min_bic);
    if min_bic > 0
        flag = 1;
    end
    if flag == 0
        % merge states min_i and min_j
        disp('merging two states...');
        if min_i > min_j
            min_s = min_j;
            max_s = min_i;
        else
            min_s = min_i;
            max_s = min_j;
        end
        
        % Model both the states using same gaussian
        count = 1;
        clear tempData;
        for k = 1:featureCount
            if stateSequence(k, 1) == min_i || stateSequence(k, 1) == min_j
                tempData(count, :) = data(k, :);
                count = count + 1;
            end
        end
        mu = mean(tempData);
        sigma = cov(tempData);
        
        % change state sequence
        for i = 1:featureCount
            if stateSequence(i, 1) == max_s
                stateSequence(i, 1) = min_s;
            end
        end
        
        for i = 1:featureCount
            if stateSequence(i, 1) > max_s
                stateSequence(i, 1) = stateSequence(i, 1) - 1;
            end
        end
        % copy new mean and covariacne matrix into allmeans and allvars        
        allStateMeans(min_s, :) = mu;
        allStateVars(min_s, :, :) = sigma;
        numElemEachState(min_s, 1) = numElemEachState(min_s, 1) + numElemEachState(max_s, 1);
        % decrease number of states by one
        numStates = numStates  - 1;
        
        % now copy all means and vars
        for i = max_s:numStates
            allStateMeans(i, :) = allStateMeans(i+1, :);
            allStateVars(i, :, :) = allStateVars(i+1, :, :);
            numElemEachState(i, 1) = numElemEachState(i+1, 1);
        end
                
        % calculate log_probability for each feature vector 
        for s = 1: numStates
            sigma = squeeze(allStateVars(s, :, :));
            prob = mvnpdf(voicedData, allStateMeans(s, :), sigma);
            prob = log(prob);
            log_prob(s, :) = prob;            
        end
        % calculate pi 
        for s = 1:numStates
            pi(s, 1) = log(numElemEachState(s, 1) / featureCount);
        end        
    else
        disp('delta bic is greater than zero, we need to stop now ');
        iterFlag = 0;        
    end
    mergeIter = mergeIter + 1;
end

% write RTTM file 
%==========================================================================
% writing the final data in RTTM file in the format specified by NIST
% allData , allFeatureFlag, numStates, fileName
file_name = 'IS1007a';
extension = '.rttm';
rttm_file = fopen(strcat(file_name, extension), 'w');

allStateSequence = zeros(allFeatureCount, 1);
for k = 1:allFeatureCount(1,1)
    if allFeatureFlag(k, 2) == 0
        allStateSequence(k, 1) = 0;
    else
        index = allFeatureFlag(k, 2);
        state = stateSequence(index, 1);
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

    

