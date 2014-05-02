%MOGFINDBESTPAIR Find the best pair of Mixture of Gaussian models
%
%     [COMB,LL,MODELS] = MOGFINDBESTPAIR(MOGS,X,I,STOPCRIT,REG)
%
% Find the two most similar Mixture of Gaussian models in the cell-array
% MOGS of Mixture of Gaussian models. The indices of the two most
% similar models are found by first computing all log-likelihoods of the
% models on their data (the data of model 'n' is given by X((I==n),:)).
% Next, all models are pairwise combined, and retrained on the combined
% data (that is, X((I==n1),:) combined with X((I==n2),:) ). The log-
% likelihood of the combined model is compared with the sum of the
% loglikelihoods of the original models. The pair for which the
% loglikelihood improvement is the largest, is the pair that is returned
% in COMB. The corresponding loglikelihood improvement is returned in
% LL. In MODELS all pairwise models are returned.
%
% See also: mogcat, moginit, hmmcombinestates

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [comb,mx,models] = mogfindbestpair(mogs,x,pth,stopcrit,reg)

% compute the ll of all combinations
N = length(mogs);
models = cell(N,N);
ll = zeros(N,N);
% first the diagonal elements:
for i=1:N
	I = find(pth==i);
	models{i,i} = mogs{i};
	ll(i,i) = sum(moglogp(x(I,:),models{i,i}));
end
% now the off-diagonals:
for i=1:N
	I = find(pth==i);
	for j=i+1:N
		J = find(pth==j);
		models{i,j} = mogcat(models{i,i},models{j,j});
		models{i,j} = mogem(x([I;J],:),models{i,j},stopcrit,reg);
		ll(i,j) = sum(moglogp(x([I;J],:),models{i,j}));
	end
end

% find the best combination:
diagll = diag(ll);
df = ll - repmat(diagll,1,N) - repmat(diagll',N,1);
% should I set the diagnoal to -inf?
df(1:(N+1):end) = -inf;
[mx,mi] = max(df,[],1);
[mx,mj] = max(mx,[],2);
mi = mi(mj);
% in the right order:
if mj<mi, tmpm=mi; mi=mj; mj=tmpm; end
comb = [mi mj];

return
