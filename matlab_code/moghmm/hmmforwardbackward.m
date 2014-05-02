%
% [GAMMA,GAMMA1,SUMXI,LL] = HMMFORWARDBACKWARD(X,MODEL,P,P_K)
%
% One forward-backward pass through a HMM model with a Mixture of
% Gaussians in each of the states.
% This is a support function for HMMEM.
%
% See also: hmmlogforwardbackward, hmmem, mdforwardbackward

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [gamma,gamma1,sumxi,ll,g] = hmmforwardbackward(x,model,p,p_k)
if nargin<3
	[p,p_k] = hmmp(x,model);
end

%setup parameters
Q = length(model.prior);
T = size(x,1);
alpha = zeros(T,Q);
beta = zeros(T,Q);
g = zeros(T,Q);
m = zeros(Q,1);
gamma = cell(Q,1);
for i=1:Q
	m(i) = length(model.pdf{i}.prior);
	gamma{i} = zeros(T,m(i));
	pi = p(:,i); pi(pi==0) = 1;
	p_k{i} = p_k{i}./repmat(pi,1,m(i));
end
scale = zeros(T,1);
sumxi = zeros(Q,Q); %makes life simpler later on

% forward initialization:
alpha(1,:) = model.prior'.*p(1,:);
[alpha(1,:),scale(1)] = vecunitnorm(alpha(1,:));
% forward induction:
for i=2:T
	alpha(i,:) = alpha(i-1,:)*model.trans.*p(i,:);
	[alpha(i,:),scale(i)] = vecunitnorm(alpha(i,:));
end
% forward termination:
if any(scale==0)
	ll = -inf;
else
	ll = sum(log(scale));
end

% backward initialization
beta(T,:) = 1;
g(T,:) = vecunitnorm(alpha(T,:).*beta(T,:));
for j=1:Q
	gamma{j}(T,:) = p_k{j}(T,:) * g(T,j);
end
% backward induction
for i=(T-1):-1:1
	%beta(i,:) = (beta(i+1,:).*p(i+1,:))*model.trans';
	bpi = (beta(i+1,:).*p(i+1,:));
	beta(i,:) = bpi*model.trans';
	beta(i,:) = vecunitnorm(beta(i,:));
	g(i,:) = vecunitnorm(alpha(i,:).*beta(i,:));
	for j=1:Q
		gamma{j}(i,:) = p_k{j}(i,:) * g(i,j);
	end
	% for the EM algorithm we need the sum over xi over all t
	sumxi = sumxi + unitnorm(model.trans.*(alpha(i,:)'*bpi));
end
gamma1 = g(1,:);

