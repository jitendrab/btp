%HMMLOGEM EM algorithm on Hidden Markov Model
%
%     MODEL = HMMLOGEM(X,MODEL,STOPCRIT,REG)
%
% Estimate the parameters given the HMM model MODEL and the data X using
% the EM algorithm. This is repeated for maximally STOPCRIT.MAXITER
% steps. The covariance matrices inside the HMM model are regularized by
% adding REG to the diagonal. This implementation uses logarithms for
% all the probabilities (numerically more stable, but also a bit slower)
%
% See also: hmmp, mogp, hmmem

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function model = hmmlogem(x,model,stopcrit,reg)
if nargin<4
	reg = 0.001;
end
if nargin<3
	stopcrit.maxiters = 10;
	stopcrit.minllimpr = 1e-5;
end

% first check: do we have a minimum duration constraint?
if isfield(model,'md') %& (model.md>1)
	usemd = 1;
	Q = length(model.prior)/model.md;
else
	usemd = 0;
	Q = length(model.prior);
end

% initialization:
[T,dim] = size(x);
iter = 1;
LL2 = -inf;

while (iter<=stopcrit.maxiters)

	% compute probabilities:
	[logp,logp_k] = hmmlogp(x,model);

	% forward backward:
	LL1 = LL2;
	if usemd
		[gamma,gamma1,sumxi,LL2] = mdhmmlogforwardbackward(x,model,logp,logp_k);
	else
		[gamma,gamma1,sumxi,LL2] = hmmlogforwardbackward(x,model,logp,logp_k);
	end
	%dd_message(4,'HMM Iteration %d (ll=%f)\n',iter,LL2);

	% is it going in the right direction?
	if (((LL2-LL1)/abs(LL1))<stopcrit.minllimpr) break; end

	% estimate the model parameters again:
	% prior and transtition probabilities
	model.prior = gamma1';
	model.trans = sparse(unitnorm(sumxi',2));
	% state models (prior, mean, cov.matrix)
	for i=1:Q
		Z = sum(gamma{i},1);
		Z(Z==0) = 1;
		mi = length(Z);
		model.pdf{i,1}.prior = vecunitnorm(Z)';
		for j=1:mi
			biggamma = repmat(gamma{i}(:,j),1,dim);
			model.pdf{i,1}.mean(j,:) = sum(x.*biggamma,1)/Z(j);
			dx = x - repmat(model.pdf{i}.mean(j,:),T,1); %use new means
			newcov = (dx.*biggamma)'*dx/Z(j) + reg*eye(dim); % regul. cov. matrix
			model.pdf{i,1}.icov(j,:,:) = inv(newcov);
		end
	end

	% how far are we?
	iter = iter+1;

end

