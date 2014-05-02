%MOG2PRTOOLS Convert MoG to Prtools mapping
%
%       W = MOG2PRTOOLS(MODEL)
%
% Convert a Mixture of Gaussians to a Prtools mapping.

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = mog2prtools(model)

[k,dim] = size(model.mean);
labs = ones(k,1);
dat.mean = model.mean;
for i=1:k
	dat.cov(:,:,i) = inv(squeeze(model.icov(i,:,:)));
end
dat.prior = model.prior';
dat.nlab = labs;

w = mapping('normal_map','trained',dat,labs,dim,1);


