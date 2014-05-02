%MOGP Probability of Mix.of Gaussians
%
%      [P,P_K] = MOGP(X,MOG)
%
% Estimate the probability P of data X of the mixture model MOG. It is
% assumed that data X has N objects in D dimensions, therefore the data
% matrix is NxD.
% The result is the probability for each of the objects X, stored in P,
% and the probability for each of the clusters (times cluster prior) in
% MOG for each of the objects, stored in P_K.
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
%See also: moglogp

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function [p,p_k] = mogp(x,mog)

K = size(mog.mean,1);
[n,dim] = size(x);
p_k = zeros(n,K);

for i=1:K
	xm = x - repmat(mog.mean(i,:),n,1);
	C = squeeze(mog.icov(i,:,:));
	d2 = sum((xm*C).*xm,2);
	p_k(:,i) = mog.prior(i)*exp(-d2/2)*sqrt(det(C)/(2*pi)^dim);
end
p = sum(p_k,2);


