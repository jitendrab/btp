function [path,fullpath,sortofprob,scale] = hmmviterbi(x,model,b)
%   PATH = HMMVITERBI(X,MODEL)
%
% Return the most likely sequence of states of hidden markov model
% MODEL for sequence X. 
%
%   PATH = HMMVITERBI(X,MODEL,B)
%
% When you already know the state and cluster probabilities B, you can
% already supply it to avoid too much overhead.
%
%   [PATH,MD_PATH,SORTOFPROB] = HMMVITERBI(...)
%
% The output MD_PATH gives the state numbers of a minimum duration
% constrained HMM, where not the full state sequence is given, but only
% the 'main' states.
% The output SORTOFPROB gives for each of the state outcomes the
% posterior.
%
% See also: hmmem, hmminit

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

if nargin<3
	b = hmmp(x,model);
end
% Take special care for the minimum duration HMM:
% copy the probabilities several times 
disp('reshaping B(posterior prob matrix),,,');
if isfield(model,'md')
	b = reshape( repmat(b,model.md,1), size(b,1), size(b,2)*model.md);
end


T = size(x,1);
Q = length(model.prior);
delta = zeros(T,Q);
psi = zeros(T,Q);
path = zeros(T,1);
scale = zeros(T,1);

% first step:
delta(1,:) = vecunitnorm(model.prior'.*b(1,:));
% next:
disp('calculating delta and psi matrix...');
for i=2:T
    disp(i);
	for j=1:Q % we have to run over the states to find the max:
		[delta(i,j), psi(i,j)] = max( delta(i-1,:) .* model.trans(:,j)' );
	end
	delta(i,:) = delta(i,:).*b(i,:);
	[delta(i,:),scale(i)] = vecunitnorm(delta(i,:));
end

% find the path back:
disp('finding backward path...');
[p,path(T)] = max(delta(T,:));
for i=(T-1):-1:1
    disp(i);
	path(i) = psi(i+1,path(i+1));
end

% condense the path for the minimum duration:
disp('condensing path for min-duration..');
if isfield(model,'md')
	fullpath = path;
	path = ceil(path/model.md);
else
	fullpath = [];
end

% give some sort of probability for each of the state outcomes:
if nargout>2
	%sortofprob = delta((path-1)*T + (1:T)');
	sortofprob = delta;
end

