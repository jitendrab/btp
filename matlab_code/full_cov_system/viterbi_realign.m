function [ stateSequence, numElemEachState, delta, si] = viterbi_realign( log_prob, pi, numStates, featureCount, numElemEachState, delta, si)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% perform viterbi realignment 
%[stateSequence] = viterbiRealignment(numStates, stateSequence, allStateMeans, allStateVars);

% step:0 Preprocessing
% take logarithm of Initial state probabilities(Pi), log_prob(prob of
% observing symbol ot at time t), a transition state probabilities


% pi is already in log
% log_prob is also in logarithm
A = log(1/numStates);
for t = 1:featureCount
    for s = 1:numStates
        delta(s, t) = 0;
        si(s, t) = 0;
    end
end

% step:1 INITIALIZATION
for i = 1: numStates
    delta(i, 1) = pi(i, 1) + log_prob(i, 1);
    si(i, 1) = 0;
end

% step: 2 recursion
for t=2:featureCount
    for j = 1:numStates
        %delta(j, t) = max
        max = -9999999;
        max_s = 1;
        for i = 1:numStates
            if max < delta(i, t-1) + A
                max = delta(i, t-1) + A;
                max_s = i;
            end
        end
        delta(j, t) = max + log_prob(j, t);
        %disp(delta(j, t));
        si(j, t) = max_s;
    end
end

% STEP:4 Termination
max = -9999999;
max_s = 1;
for i = 1:numStates
    if max < delta(i, featureCount)
        max = delta(i, featureCount);
        max_s = i;
    end
end
P = max;
q = max_s;
stateSequence(featureCount, 1) = q;
% step:5  backtracking
for t = featureCount-1:-1:1
    ind = stateSequence(t+1, 1);
    stateSequence(t, 1) = si(ind, t+1);
end

% change pi and numElemEachState
for i = 1:featureCount
    ind = stateSequence(i, 1);
    numElemEachState(ind, 1) = numElemEachState(ind, 1) + 1;
end

for s = 1:numStates
    pi(s, 1) = log(numElemEachState(s, 1) / featureCount);
end

end