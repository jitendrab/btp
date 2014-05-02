%   [P,P_K] = HMMP(X,MODEL)
%
% Compute the state and cluster probabilities P and P_K for dataset X on
% Hidden Markov model MODEL. Note that the transition probabilities are
% NOT taken into account.
%
% See also: hmmlogp, mogp, hmmem

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [p,p_k] = hmmp(x,model)

if isfield(model,'md') && (model.md>1)
	% more complicated: we have minimum duration constraints...
	Q = length(model.prior)/model.md;
else
	% standard HMM, without minimum duration constraints
	Q = length(model.prior);
end

p = zeros(size(x,1),Q);
p_k = cell(1,Q);
for j=1:Q
	[p(:,j), p_k{j}] = mogp(x,model.pdf{j});
end

