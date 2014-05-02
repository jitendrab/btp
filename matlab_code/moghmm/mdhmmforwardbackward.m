% [gamma,gamma1,sumxi,ll] = mdhmmforwardbackward(x,model,p,p_k)
%
% One forward-backward pass through a minimum-duration HMM model with a
% Mixture of Gaussians in each of the states.
% This is a support function for HMMEM.
%
% See also: hmmforwardbackward, hmmem

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [gamma,gamma1,sumxi,ll,g] = mdhmmforwardbackward(x,model,p,p_k)
if nargin<3
	[p,p_k] = hmmp(x,model);
end

%setup parameters
Qmd = length(model.prior);
Q = Qmd/model.md;
T = size(x,1);
alpha = zeros(T,Qmd);
beta = zeros(T,Qmd);
g = zeros(T,Qmd);
m = zeros(Q,1);
gamma = cell(Q,1);
for i=1:Q
	m(i) = length(model.pdf{i}.prior);
	gamma{i} = zeros(T,m(i));
	pi = p(:,i); pi(pi==0) = 1;
	p_k{i} = p_k{i}./repmat(pi,1,m(i));
end
scale = zeros(T,1);
sumxi = zeros(Qmd,Qmd); %makes life simpler later on
% expand the probabilities p, to make multiplication later easier
% (although it requires now quite some memory...)
p = reshape( repmat(p,model.md,1), size(p,1), size(p,2)*model.md);

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
	Jsubstates = ((j-1)*model.md+1):(j*model.md);
	gamma{j}(T,:) = p_k{j}(T,:) * sum(g(T,Jsubstates),2);
end
% backward induction
for i=(T-1):-1:1
	beta(i,:) = (beta(i+1,:).*p(i+1,:))*model.trans';
	beta(i,:) = vecunitnorm(beta(i,:));
	g(i,:) = vecunitnorm(alpha(i,:).*beta(i,:));
	for j=1:Q
		Jsubstates = ((j-1)*model.md+1):(j*model.md);
		gamma{j}(i,:) = p_k{j}(i,:) * sum(g(i,Jsubstates),2);
	end
	% for the EM algorithm we need the sum over xi over all t
	sumxi = sumxi + ...
		unitnorm(model.trans'.*(alpha(i,:)'*(beta(i+1).*p(i+1,:))));
end
gamma1 = g(1,:);

