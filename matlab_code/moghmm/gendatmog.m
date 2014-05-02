%GENDATMOG Generate data from a MOG 
%
%     X = GENDATMOG(N,MOG)
%
% Generate N objects X from a Mixture of Gaussians model MOG.
% The model should contain the following fields:
%    MOG.prior    Kx1 cluster priors
%    MOG.mean     KxD cluster means
%    MOG.icov     KxDxD inverse covariance matrices
%
% See also: gendatmoghmm

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function x = gendatmog(n,mog)

% Initialize:
N = rand(n,1);
x = zeros(n,size(mog.mean,2));
K = length(mog.prior);

% First decide from which mixture component each object is generated:
cump = cumsum(mog.prior);
I{1} = find(N<cump(1));
for i=2:K
	I{i} = find((N>=cump(i-1))&(N<cump(i)));
end

% Next generate data from the components:
for i=1:K
	if length(I{i})>0 %generate data when more than 0 objects are req.
		sig = inv(squeeze(mog.icov(i,:,:)));
		x(I{i},:) = randmvn(length(I{i}), mog.mean(i,:), sig);
	end
end

