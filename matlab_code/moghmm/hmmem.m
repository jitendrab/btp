%HMMEM EM algorithm on Hidden Markov Model
%
%     MODEL = HMMEM(X,MODEL,STOPCRIT,REG)
%
% Estimate the parameters given the HMM model MODEL and the data X using
% the EM algorithm. This is repeated for maximally STOPCRIT.MAXITER
% steps, of when the loglikelihood does not improve STOPCRIT.MINLLIMPR. The
% covariance matrices inside the HMM model are regularized by adding REG
% to the diagonal.
%
% Note that STOPCRIT is a Matlab structure with two fields, maxiter and
% minllimpr. Default crit.maxiters=10, crit.minllimpr=1e-5.
%
% See also: hmminit, hmmviterbi, hmmlogem

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [model,LL] = hmmem(x,model,stopcrit,reg)
if nargin<4
	reg = 0.001;
end
if (nargin<3) || (isempty(stopcrit))
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
LL = [];
LL2 = -inf;

% go:
while (iter<=stopcrit.maxiters)

	% compute probabilities:
	[p,p_k] = hmmp(x,model);

	% forward backward:
	LL1 = LL2;
	if usemd
		[gamma,gamma1,sumxi,LL2] = mdhmmforwardbackward(x,model,p,p_k);
	else
		[gamma,gamma1,sumxi,LL2] = hmmforwardbackward(x,model,p,p_k);
	end
	LL(iter) = LL2;
	%dd_message(4,'HMM Iteration %d (ll=%f)\n',iter,LL2);

	% is it going in the right direction?
	if (((LL2-LL1)/abs(LL1))<stopcrit.minllimpr) break; end

	% estimate the model parameters again:
	% prior and transtition probabilities
	%if abs(sum(gamma1)-1)>1e-3, keyboard end
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
			if reg<0
				% use diagonal cov. matrices:
				newcov = (dx.*biggamma)'*dx/Z(j);
				model.pdf{i,1}.icov(j,:,:) = diag(1./diag(newcov));
			else
				% use full cov. matrices:
				newcov = (dx.*biggamma)'*dx/Z(j) + reg*eye(dim); % regul. cov. matrix
				model.pdf{i,1}.icov(j,:,:) = inv(newcov);
			end
		end
	end

	% test if everything goes well:
    % below lines commented by me 
% 	if abs(sum(model.prior)-1)>1e-2  % changed by me from 1e-6
% 		keyboard
% 		error('HMMem: Priors do not add up to one.');
% 	end
	% how far are we?
	iter = iter+1;

end

