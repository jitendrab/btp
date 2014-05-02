%MOGEM Apply EM on Mixture of Gaussian
%
%      MOG = MOGEM(X,MOG,STOPCRIT,REG)
%
% Estimate a Mixture of Gaussians model MOG on dataset X using the EM
% algorithm. The EM algorithm is run for STOPCRIT.MAXITER times or when
% the likelihood does not improve more than STOPCRIT.MINLLIMPR. Each of
% the covariance matrices is regularized by adding REG to the diagonal.
%
% Default we use:  STOPCRIT.maxiters = 10
%                  STOPCRIT.minllimpr = 1e-5
%                  REG = 1e-5;
% See also: mogp, moginit

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function mog = mogem(x,mog,stopcrit,reg)

if nargin<4
	reg = 1e-5;
end
if nargin<3
	stopcrit.maxiters = 10;
	stopcrit.minllimpr = 1e-5;
end
[n,dim] = size(x);
k = length(mog.prior);
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
	new_priors = sum(p_k,1)';
	new_priors(new_priors==0) = 10*realmin;
	new_mean = p_k' * x;

	% and make a new model
	mog.prior = new_priors/n;
	mog.mean = bsxfun(@rdivide, new_mean, new_priors);
	for i=1:k
		df = bsxfun(@times, bsxfun(@minus,x,mog.mean(i,:)),sqrt(p_k(:,i)));
		mog.icov(i,:,:) = inv((df'*df)/new_priors(i) + reg*eye(dim));
	end

	% how far are we?
	iter = iter+1;
end

