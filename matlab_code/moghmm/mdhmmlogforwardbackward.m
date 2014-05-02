function [gamma,gamma1,sumxi,ll,logg] = mdhmmlogforwardbackward(x,model,logp,logp_k)
% [gamma,gamma1,sumxi,ll] = mdhmmlogforwardbackward(x,model,logp,logp_k)
%
% One forward-backward pass through a minimum-duration HMM model with a
% Mixture of Gaussians in each of the states.
% This version uses the logarithm of all the probabilities, to ensure a
% bit better numerical behaviour.
if nargin<3
	[logp,logp_k] = hmmlogp(x,model);
end

% Setup parameters
% Note that the total number of states is now much larger: 
% instead of Q (= number of PDFs) we have to extend each state by MD
% substates. Therefore the total number of states is Q*MD
Qmd = length(model.prior);
Q = Qmd/model.md;
T = size(x,1);
% Allocate:
warning off MATLAB:log:logOfZero;
	logA = log(model.trans);
warning on MATLAB:log:logOfZero;
logalpha = zeros(T,Qmd);
logbeta = zeros(T,Qmd);
logg = zeros(T,Qmd);
m = zeros(Q,1);
gamma = cell(Q,1);
for i=1:Q
	m(i) = length(model.pdf{i}.prior);
	gamma{i} = zeros(T,m(i));
	logpi = logp(:,i);
	logp_k{i} = logp_k{i} - repmat(logpi,1,m(i));
end
sumxi = zeros(Qmd,Qmd); %makes life simpler later on
% We assume that the substates of a state all have the save PDF.
% Therefore the probabilities p can be copied to each substate. We have
% MD substates per state, so we expand the probabilities p from each
% state and copy it to each substate
% (although it requires now quite some memory...)
logp = reshape( repmat(logp,model.md,1), size(logp,1), size(logp,2)*model.md);

% forward initialization:
warning off MATLAB:log:logOfZero;
	logalpha(1,:) = log(model.prior)' + logp(1,:);
warning on MATLAB:log:logOfZero;
% forward induction:
for i=2:T
	%for j=1:Qmd
	%	logalpha(i,j) = hmmlogsum(logalpha(i-1,:)+logA(:,j)',2) + logp(i,j);
	%end
	logalpha(i,:) = hmmlogsum(repmat(logalpha(i-1,:)',1,Qmd)+logA,1)+logp(i,:);
end
% forward termination:
ll = sum(logalpha(end,:));

% backward initialization
logg(T,:) = logalpha(T,:) - hmmlogsum(logalpha(T,:),2);
for j=1:Q
	% find which substates share a pdf:
	Jsubstates = ((j-1)*model.md+1):(j*model.md);
	% Originally, we had a gamma per state. Now the states consists of
	% substates, and we sum the substate probabilities
	gamma{j}(T,:) = exp(logp_k{j}(T,:) + hmmlogsum(logg(T,Jsubstates),2));
end

% backward induction
for i=(T-1):-1:1
	for j=1:Qmd
		logbeta(i,j) = hmmlogsum(logA(j,:)+logbeta(i+1,:)+logp(i+1,:),2);
	end
	logg(i,:) = logalpha(i,:) + logbeta(i,:);
	logg(i,:) = logg(i,:) - hmmlogsum(logg(i,:),2);
	for j=1:Q
		Jsubstates = ((j-1)*model.md+1):(j*model.md);
		gamma{j}(i,:) = exp(logp_k{j}(i,:) + hmmlogsum(logg(i,Jsubstates),2));
	end
	% for the EM algorithm we need the sum over xi over all t
	lognewxi = repmat(logalpha(i,:)',1,Qmd) + logA + ...
		repmat(logbeta(i+1,:)+logp(i+1,:),Qmd,1);
	sumxi = sumxi + exp(lognewxi - hmmlogsum(lognewxi(:),1));
end

% annoying numerics:
for j=1:Q
	gamma{j}(isnan(gamma{j})) = 0.0;
end
gamma1 = sum(exp(lognewxi - hmmlogsum(lognewxi(:))),2)';

