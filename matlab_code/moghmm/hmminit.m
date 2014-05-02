%HMMINIT Initialize HMM model
%
%      MODEL = HMMINIT(N,X)
%      MODEL = HMMINIT(N,DIM,MSPREAD,VARS)
%
% Create a Hidden Markov Model where each of the states is a Mixture of
% Gaussians. The vector N contains the number of clusters in each of the
% states (therefore length(n) = # states). When a dataset X is supplied,
% the dimensionality, the spread in the means and the initial covariance
% matrices are derived from this dataset X. Otherwise, DIM contains the
% dimensionality of the observables, MSPREAD gives the variance of the
% Gaussian from which the initial means are drawn, and VARS gives the
% variance of the initial covariance matrices. The state priors and the state
% transitions are uniform.
%
% So, W=initmoghmm([4 1 5],2) generates a HMM with 3 states, where the
% first state is modeled by a 4-cluster MoG in 2D. The second state is
% modeled by just a single 2D gaussian, and the third state is a MoG
% with 5 clusters.
%
% See also: hmminitwithsegments, hmmem, hmmviterbi

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function model = hmminit(n,dim,mspread,varS)
if (nargin<3) && (size(dim,1)>1)
	% a dataset X is supplied, use the cov.matrix for the generation of
	% means and cov.matrices:
	x = dim;
	xmean = mean(x);
	xcov = cov(x);
	varS = xcov;
else
	% use the user-supplied parameters:
	if nargin<4
		varS = 4;
	end
	if nargin<3
		mspread = 5;
	end
	xmean = zeros(1,dim);
	xcov = mspread*eye(dim);
	varS = varS*eye(dim);
end

% for now:
%randn('state',0);

if length(n)==1
	warning('moghmm:oneGaussianPerState', ...
        'Assuming that all states have a single Gaussian.');
	n = ones(n,1);
end
Q = length(n);
iC = inv(varS);

model.prior = repmat(1/Q,Q,1);
model.trans = repmat(1/Q,Q,Q);
for i=1:Q
	model.pdf{i,1}.prior = repmat(1/n(i),n(i),1);
	model.pdf{i,1}.mean = randmvn(n(i),xmean,xcov);
	for j=1:n(i)
		model.pdf{i,1}.icov(j,:,:) = iC;
	end
end


