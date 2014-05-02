%MOGEM_WEIGHTED Apply EM on Mixture of Gaussian on weighted data
%
%      MOG = MOGEM_WEIGHTED(X,W,MOG,STOPCRIT,REG)
%
% Estimate a Mixture of Gaussians model MOG on dataset X using the EM
% algorithm. In W the weights for each of the objects in X are given.
% The EM algorithm is run for STOPCRIT.MAXITER times or when the
% likelihood does not improve more than STOPCRIT.MINLLIMPR. Each of the
% covariance matrices is regularized by adding REG to the diagonal.
%
% Default we use:  STOPCRIT.maxiters = 10
%                  STOPCRIT.minllimpr = 1e-5
%                  REG = 1e-5;
% See also: mogp, moginit, mogem

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function mog = mogem_weighted(x,w,mog,stopcrit,reg)

if nargin<5
	reg = 1e-5;
end
if nargin<4
	stopcrit.maxiters = 10;
	stopcrit.minllimpr = 1e-5;
end
[n,dim] = size(x);
k = length(mog.prior);
if size(w,1)~=size(x)
	error('Number of weights should be equal to the number of objects.');
end
if (any(w<0))
	error('Please use positive weights.');
end
w = w/sum(w); % normalize the weights

iter = 1;
LL2 = -inf;

while (iter<=stopcrit.maxiters)

	% estimate responsibilities
	[p,p_k] = mogp(x,mog);

	% do we improve?
	LL1 = LL2;
	LL2 = sum(log(p));
	if ((LL2-LL1)/abs(LL1))<stopcrit.minllimpr break; end

	% normalize:
	p(p==0) = 1; % don't divide by zero!
	p_k = bsxfun(@rdivide,p_k,p);

	% update
	new_priors = sum( bsxfun(@times,p_k,w), 1)';
	new_priors(new_priors==0) = 10*realmin;
	mog.prior = new_priors/n;

	% and make a new model
	new_mean = bsxfun(@times, p_k, w)' * x;
	mog.mean = bsxfun(@rdivide,new_mean,new_priors);
	for i=1:k
		df = bsxfun(@times, bsxfun(@minus,x,mog.mean(i,:)), sqrt(w.*p_k(:,i)));
		mog.icov(i,:,:) = inv((df'*df)/new_priors(i) + reg*eye(dim));
	end

	% how far are we?
	iter = iter+1;
end

