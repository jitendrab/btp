%
% [GAMMA,GAMMA1,SUMXI,LL] = HMMLOGFORWARDBACKWARD(X,MODEL,P,P_K)
%
% One forward-backward pass through a HMM model with a Mixture of
% Gaussians in each of the states.
% This is a support function for HMMEM.
%
% See also: hmmforwardbackward, hmmem, mdhmmforwardbackward

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [gamma,gamma1,sumxi,ll,logg] = hmmlogforwardbackward(x,model,logp,logp_k)
if nargin<3
	[logp,logp_k] = hmmlogp(x,model);
end

%setup parameters
warning off MATLAB:log:logOfZero;
	logA = log(model.trans);
warning on MATLAB:log:logOfZero;
Q = length(model.prior);
T = size(x,1);
logalpha = zeros(T,Q);
logbeta = zeros(T,Q);
logg = zeros(T,Q);

% also do this here:
m = zeros(Q,1);
gamma = cell(Q,1);
for i=1:Q
	m(i) = length(model.pdf{i}.prior);
	gamma{i} = zeros(T,m(i));
	% normalize the mixture pdf in p_k
	logpi = logp(:,i);
	logp_k{i} = logp_k{i} - repmat(logpi,1,m(i));
end
sumxi = zeros(Q,Q); %makes life simpler later on

% forward initialization:
warning off MATLAB:log:logOfZero;
	logalpha(1,:) = log(model.prior)' + logp(1,:);
warning on MATLAB:log:logOfZero;
% forward induction:
for i=2:T
	%for j=1:Q
	%	logalpha(i,j) = hmmlogsum(logalpha(i-1,:)+logA(:,j)',2) + logp(i,j);
	%end
	logalpha(i,:) = hmmlogsum(repmat(logalpha(i-1,:)',1,Q)+logA,1)+logp(i,:);
end
% gives the likelihood:
ll = sum(logalpha(end,:));

% backward initialization
logg(T,:) = logalpha(end,:) - hmmlogsum(logalpha(end,:),2);
for j=1:Q
	gamma{j}(T,:) = exp(logp_k{j}(T,:) + logg(T,j));
end
% backward induction
for i=(T-1):-1:1
	for j=1:Q
		logbeta(i,j) = hmmlogsum(logA(j,:)+logbeta(i+1,:)+logp(i+1,:),2);
	end
	% computation of log(gamma) is now possible (called logg here):
	logg(i,:) = logalpha(i,:) + logbeta(i,:);
	logg(i,:) = logg(i,:) - hmmlogsum(logg(i,:),2);
	% finally, the gamma_k is computed (called gamma here):
	for j=1:Q
		gamma{j}(i,:) = exp(logp_k{j}(i,:) + logg(i,j));
	end
	% for the EM algorithm we need the sum over xi over all t
	%keyboard
	lognewxi = repmat(logalpha(i,:)',1,Q) + logA + ...
		repmat(logbeta(i+1,:)+logp(i+1,:),Q,1);

	sumxi = sumxi + exp(lognewxi - hmmlogsum(lognewxi(:)));
end

% annoying numerics:
for j=1:Q
	gamma{j}(isnan(gamma{j})) = 0.0;
end
gamma1 = sum(exp(lognewxi - hmmlogsum(lognewxi(:))),2)';

