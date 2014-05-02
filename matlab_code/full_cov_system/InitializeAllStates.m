function [allStateMeans, allStateVars] = InitializeAllStates(data, INITIAL_NUM_STATES, allStateMeans, allStateVars, featureCount)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
j = 1;
block = featureCount/ INITIAL_NUM_STATES;
x = zeros(block, 19);
w = zeros(block, 1);
for s=1:INITIAL_NUM_STATES
    i = 1;
    for j=(s-1)*block+1:s*block
        x(i, :) = data(j, :);
        w(i, 1) = 1;
        i = i+1;
    end
    disp(size(x));
    [mean, variance, weight, avg_prob, fisher, log_prob, gg] = gaussmix(x, [], [], 1, 'v', w, w);
    % we should store obtained mean and variance into all mixture means and
    % variances
    allStateMeans(s, :) = mean(1, :);
    allStateVars(s, :, :) = variance(:, :);
end
end
