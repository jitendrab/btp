%HMMINITWITHSEGMENTS Initialize HMM model
%
%      HMM = HMMINITWITHSEGMENTS(X,N)
%
% Split the timeseries X into N segments and fit one state with one
% Gaussian on each of the segments.
%
% See also: hmminit, moginit

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function model = hmminitwithsegments(x,Q)

[n,dim] = size(x);
% the simple things in the model:
model.prior = repmat(1/Q,Q,1);
model.trans = repmat(1/Q,Q,Q);
% split it, and train a hidden node on each of the segments:
I = ceil((1:n)*Q/n);
for i=1:Q
	ai = x(find(I==i),:);
	model.pdf{i}.prior = 1;
	model.pdf{i}.mean = mean(ai);
	c = inv(cov(ai));
	model.pdf{i}.icov = reshape(c,[1,dim,dim]);
end



