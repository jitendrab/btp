%HMMCOMBINECLUSTERS
%
%   [MOGS,I,DELTALL] = HMMCOMBINECLUSTERS(MOGS,X,STOPCRIT,REG)
%
% Combine two Mixture of Gaussian models from the cellarray MOGS into
% one by comparing the Mixture of Gaussian probability models for each
% of the clusters, and checking if the combination of the two clusters
% gives a higher log-likelihood on the data X.
function [mogs,I,mx] = hmmcombineclusters(mogs,x,stopcrit,reg)

% set parameters
[m,dim] = size(x);
N = length(mogs);

while (N>1)
	% find out to which MoG each object is assigned:
	p = zeros(m,N);
	for i=1:N
		p(:,i) = moglogp(x,mogs{i});
	end
	[dummy,I] = max(p,[],2);

	% remove the models for which none of the objects is assigned:
	for i=N:-1:1
		if isempty(find(I==i))
			dd_message(2,'Remove model %d\n',i);
			mogs{i} = [];
			J = find(I>i);
			I(J) = I(J)-1;
		end
	end

	% find the best combination of MOGs
	[comb,mx,models] = mogfindbestpair(mogs,x,I,stopcrit,reg);
	if (comb(1)==comb(2)) || (mx<0)
		dd_message(3,'To combine [%d,%d] is not possible (deltall=%f).\n',...
			comb(1),comb(2),mx);
		mx = -abs(mx); %make sure it is negative
		break;
	else
		dd_message(3,'Combine [%d,%d]->%d (deltall=%f).\n',comb(1),comb(2),comb(1),mx);
	end

	% keep the combined model, remove the old model:
	N = size(mogs,1);
	mogs = mogs([1:N+1:end]);
	%mogs = diag(models);
	mogs{comb(1)} = models{comb(1),comb(2)};
	mogs(comb(2)) = [];
	N = length(mogs);
end


