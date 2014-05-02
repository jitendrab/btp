%MOGLOGP Log-Probability of Mix.of Gaussians
%
%      [LOGP,LOGP_K] = MOGLOGP(X,MOG)
%
% Estimate the log-probability P of data X of the mixture model MOG. It
% is assumed that data X has N objects in D dimensions, therefore the
% data matrix is NxD.
% The result is the log-probability for each of the objects X, stored in
% LOGP, and the log-probability for each of the clusters in MOG for each
% of the objects, LOGP_K.
%
% The mixture of Gaussians is encoded in a Matlab structure. When you
% want to have K clusters in D dimensions, you have to supply:
%   mog.prior   cluster prior,                         size Kx1
%   mog.mean    cluster means,                         size KxD
%   mog.icov    (inverse) cluster covariance matrices, size KxDxD
%
% Example:    mog.prior = [1/3; 1/3; 1/3];
%             mog.mean = [0 0; 1 1; 4 0];
%             mog.icov(1,:,:) = reshape( inv(eye(2)), [1 2 2]);
%             mog.icov(2,:,:) = reshape( inv(eye(2)), [1 2 2]);
%             mog.icov(3,:,:) = reshape( inv(eye(2)), [1 2 2]);
%
% See also: mogp, moginit

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [logp,logp_k] = moglogp(x,mog)

K = size(mog.mean,1);
[n,dim] = size(x);
logp_k = zeros(n,K);

for i=1:K
	xm = x - repmat(mog.mean(i,:),n,1);
	C = squeeze(mog.icov(i,:,:));
	d2 = sum((xm*C).*xm,2);
	logp_k(:,i) = log(mog.prior(i)) -d2/2 + ...
		log(det(C))/2 - dim*log(2*pi)/2;
end
logp = hmmlogsum(logp_k,2);


