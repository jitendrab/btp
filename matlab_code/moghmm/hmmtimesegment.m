%HMMTIMESEGMENT Timeseries segmentation by HMM
%
%     [HMM,PTH,COMB] = HMMTIMESEGMENT(X,N,MD,STOPCRIT,REG)
%
% Segment the timeseries X into segments with minimum duration MD.
% Each timepoint in the timeseries (of length T) can be a feature vector
% of length DIM. The total size of X should therefore be TxDIM.
% Initially the timeseries is split into N consequitive segments, each
% modeled by a single state in an HMM using a single Gaussian pdf. The
% Gaussian models per state are merged until the total (log-)likelihood
% does not improve anymore.
%
% See also: hmmcombinestates, hmmem

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [hmm,pth,comb] = hmmtimesegment(x,N,MD,stopcrit,reg)

%dd_message(4,'Initialize HMM and train\n');
disp('Initializing HMM and training....');
% split it, and train a hidden node on each of the segments:
hmm = hmminitwithsegments(x,N);
% store it in an HMM:
hmm = hmm2mdhmm(hmm,MD);  % minimum state duration
%dd_message(4,'Find Viterbi path\n');
disp('calculating viterbi path....');
pth = hmmviterbi(x,hmm);

% now combine states till it becomes worse
disp('combining states till log likelihood increases....');
deltall=inf;
run = 0;
comb = [];
while (deltall>0)
	% and try if combining helps:
	run = run+1;
	%dd_message(2,'Run %d\n',run);
	[hmm,pth,comb(run,:),deltall] = hmmcombinestates(hmm,x,pth,stopcrit,reg);
end
%dd_message(3,'Experiment done.\n');
disp('experiment done....in timesegment function');
end

