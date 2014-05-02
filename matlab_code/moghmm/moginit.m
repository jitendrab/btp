%MOGINIT Initialize a Mixture of Gaussians model
%
%    MOG = MOGINIT(X,K)
%
% Initialize a Mixture of Gaussians with K mixture components on dataset
% X.
% The K priors are uniform 1/K.
% The K means are randomly drawn from dataset X.
% The K cov. matrices are all diagonal matrices with the same diagonal
% as the original data X.
%
% See also: hmminit, mogem

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function mog = moginit(x,k)

xvar = diag(nancov(x));
iC = diag(1./xvar);
mog.prior = repmat(1/k,k,1);
I = randperm(size(x,1));
mog.mean = x(I(1:k),:);
for i=1:k
	mog.icov(i,:,:) = iC;
end
