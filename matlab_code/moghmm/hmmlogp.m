%   [LOGP,LOGP_K] = HMMLOGP(X,MODEL)
%
% Compute the state and cluster log-probabilities LOGP and LOGP_K for
% dataset X on Hidden Markov model MODEL. Note that the transition
% probabilities are NOT taken into account.
%
% See also: hmmp, mogp, hmmem

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [logp,logp_k] = hmmlogp(x,model)

if isfield(model,'md') && (model.md>1)
	% more complicated: we have minimum duration constraints...
	Q = length(model.prior)/model.md;
else
	% standard HMM, without minimum duration constraints
	Q = length(model.prior);
end

logp = zeros(size(x,1),Q);
logp_k = cell(1,Q);
for j=1:Q
	[logp(:,j), logp_k{j}] = moglogp(x,model.pdf{j});
end

