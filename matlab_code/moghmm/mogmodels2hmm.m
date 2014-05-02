%MOGMODELS2HMM Combine mixture models to HMM model
%
%    HMM = MOGMODELS2HMM(MOG_MODELS)
%
% Combine the Mixture of Gaussian models, that are stored in the
% cell-array MOG_MODELS, into one HMM model. The state-transition matrix
% and the state priors are assumed to be uniform.
%
% See also: moginit, hmminit

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function hmm = mogmodels2hmm(models)

N = length(models);
hmm.prior = repmat(1/N,N,1);
hmm.trans = repmat(1/N,N,N);
for i=1:N
	hmm.pdf{i,1} = models{i};
end

