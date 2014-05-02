%HMM_MIL Train a HMM classifier
%
%   W = HMM_MIL(A,FRAC,N,REG,MAXITER)
%
% INPUT
%    A      MIL dataset
%    FRAC   Fraction of instances taken into account in evaluation
%    N      Vector containing the number of clusters in each of the states
%
% OUTPUT
%    W      Trained simple mapping
%
% DESCRIPTION
% Train a HMM classifier with in each of the states a Mixture of
% Gaussians. The number of Gaussians in each state is defined by N=[n_1,
% n_2, ..., N_m], where m is the total number of states.
% The covariances are regularized by adding REG to the diagonal.
% The classifier outputs for each time point the probability for each of
% the states.
%
% SEE ALSO
%   LABELBAGP

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function W = hmm_mil(a,frac,N,reg,maxiter)

% Take care of empty/not-defined arguments:
if nargin < 5, maxiter = 100; end
if nargin < 4, reg = 1e-3; end
if nargin < 3, N = [1; 1]; end
if nargin < 2, frac = 'presence'; end
if nargin < 1 || isempty(a) 
	% When no inputs are given, we are expected to return an empty
	% mapping:
	W = mapping(mfilename,{frac,N,reg,maxiter});
	W = setname(W,'HMM with m=%d',length(N));
	return
end

if ~ismapping(frac)||isuntrained(frac)           %training
	if ~hasmilbags(a)
		error('I need a MIL dataset to define the sequences.');
	end

	% get the sequences and other parameters:
	[bags,baglab,bagid,Ibag] = getbags(a);
	[n,dim,c] = getsize(a);
	crit.maxiters=1; %each sequence is used once in an interation
	crit.minllimpr = 1e-5;
	LL2 = -inf;

	%only use the positive sequences:
	I = find(ispositive(baglab));
	bags = bags(I);

	% initialize the HMM:
	hmm = hmminit(N,cell2mat(bags));
	% train:
	iter = 0;
	while (iter<=maxiter) % do several iterations
		dd_message(4,'%d/%d ',iter,maxiter);
		LL1 = LL2;
		LL2 = 0;
		for i=1:length(bags) % train all sequences:
			[hmm,LLtmp] = hmmem(bags{i},hmm,crit,reg);
			LL2 = LL2+LLtmp;
		end
		% are we good enough?
		if (((LL2-LL1)/abs(LL1))<crit.minllimpr) break; end
		iter = iter+1;
	end
	if abs(sum(hmm.prior)-1)>1e-6
		warning('moghmm:hmm_mil:zeroPriors',...
			'HMM training failed: priors do not add up to 1.');
	end

	%and save all useful data in a structure:
	W.frac = frac;  % a fraction should *always* be defined
	W.hmm = hmm;
	W = mapping(mfilename,'trained',W,getlablist(a),size(a,2),c);
	W = setname(W,'HMM with m=%d',length(N));

else                               %testing

	% Unpack the mapping and dataset:
	W = getdata(frac);
	[bags,baglab,bagid] = getbags(a);
   % Evaluate the classifier on new sequences:
	n = size(bags,1);
	out = [];
	for i=1:n
		[pth,tmp,sortofprob,scale] = hmmviterbi(bags{i},W.hmm);
		[mx,I] = max(sortofprob,[],2);
		lab = [lab;pth];
		out = [out; sortofprob];
	end
	% store:
	W = setdat(a,out,frac);
end

return

