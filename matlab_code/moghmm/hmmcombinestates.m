%HMMCOMBINESTATES
%
%   [HMM,PTH,COMB,DELTALL] = HMMCOMBINESTATES(HMM,X,PTH,STOPCRIT,REG)
%
% Combine two states from an HMM into one by comparing the Mixture of
% Gaussian probability models for each of the states, and checking if
% the combination of the two states gives a higher log-likelihood on the
% data X.
%
% See also: hmminit, hmmem, hmminitwithsegments

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [hmm,pth,comb,mx,ll] = hmmcombinestates(hmm,x,pth,stopcrit,reg)

% set parameters
N = length(hmm.pdf);
MD = hmm.md;
% remove the states for which none of the objects is assigned:
for i=N:-1:1
	if isempty(find(pth==i))
		%dd_message(2,'Remove state %d\n',i);
		disp('removing state...');
		hmm.pdf(i) = [];
		J = find(pth>i);
		pth(J) = pth(J)-1;
	end
end

% find the best combination of MOGs
%dd_message(4,'Compute all state LLs\n');
[comb,mx,models] = mogfindbestpair(hmm.pdf,x,pth,stopcrit,reg);
if (comb(1)==comb(2)) || (mx<0)
	%dd_message(3,'To combine [%d,%d] is not possible (deltall=%f).\n',...
    disp('combining is not possible...');
		%comb(1),comb(2),mx);
	mx = -abs(mx); %make sure it is negative
	return;
else
	%dd_message(3,'Combine [%d,%d]->%d (deltall=%f).\n',comb(1),comb(2),comb(1),mx);
    disp('combining two states...');
end

% keep the combined model, remove the old model and make a new HMM
% again:
N = size(models,1);
%newmodels = diag(models);
newmodels = models(1:N+1:end);
newmodels{comb(1)} = models{comb(1),comb(2)};
newmodels(comb(2)) = [];
hmm = hmm2mdhmm(mogmodels2hmm(newmodels), MD);

% train again:
%dd_message(4,'HMM train\n');
disp('train hmm again...');
%newhmm = hmmlogem(x,hmm,stopcrit,reg);
newhmm = hmmem(x,hmm,stopcrit,reg);
if all(hmm.trans(:)==0)
	warning('moghmm:allTransitionsZero',...
        'Something went wrong: all transitions are zero!');
	keyboard
end
hmm = newhmm;

% evaluate the timeseries with this model
%dd_message(4,'Find Viterbi path\n');
disp('finding viterbi path again...');
pth = hmmviterbi(x,hmm);
end

