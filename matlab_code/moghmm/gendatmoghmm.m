%GENDATMOGHMM Generate data from HMM model
%
%    [X,LABX] = GENDATMOGHMM(HMM,NRSEGM,TMIN)
%
% Generate data X from a Hidden Markov model HMM. Allow for NRSEGM state
% changes, where each state generates at least TMIN consecutive objects.
% Output variable LABX contains the 'true' states that the data was
% generated from.
%
% See also: gendatmog, randmvn

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [x,labx] = gendatmoghmm(moghmm,nrsegm,Tmin)

if nargin<3
	Tmin = 10;
end
if nargin<2
	nrsegm = 15;
end

% how long does each segment take?
T = Tmin + round(Tmin*rand(nrsegm,1));

% randomly pick a starting segment:
Q = length(moghmm.prior);
cumpr = cumsum(moghmm.prior);
pth(1) = sum(cumpr<rand(1,1))+1;

% concatenate all generated sequences:
x = [];
labx = [];
for t=1:nrsegm
	x = [x; gendatmog(T(t), moghmm.pdf{pth(t)})];
	labx = [labx; repmat(pth(t),T(t),1)];
	% next state:
	pr = zeros(Q,1); pr(pth(t))=1; pr = moghmm.trans*pr;
	cumpr = cumsum(pr);
	pth(t+1) = sum(cumpr<rand(1,1))+1;
end

