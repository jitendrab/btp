%MOGCAT Concatenate two Mixture of Gaussians
%
%       MODEL = MOGCAT(MODEL1,MODEL2)
%       MODEL = MOGCAT(MODEL1,MODEL2,WEIGHTS)
%
% Concatenate two Mixture of Gaussian models MODEL1 and MODEL2 into a
% new Mixture of Gaussian MODEL. When the WEIGHTS are not given, assume
% that each of the models has a prior of 50%.
%
%       MODEL = MOGCAT(MODELS)
%       MODEL = MOGCAT(MODELS,WEIGHTS)
%
% Concatenate the Mixture of Gaussian models stored in the cell-array
% MODELS into one Mixture of Gaussians MODEL.
%
% See also: moginit, hmmcombinestates

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function model = mogcat(model1,model2,weights)

if isa(model1,'cell')
	n = length(model1);
	if nargin<2
		weights = repmat(1/n,n,1);
	else
		weights = model2;
	end
	model = model1{1};
	model.prior = weights(1)*model.prior;
	for i=2:n
		model.prior = [model.prior(:); weights(i)*model1{i}.prior(:)];
		model.mean = [model.mean; model1{i}.mean];
		model.icov = cat(1,model.icov,model1{i}.icov);
	end
else
	if nargin<3
		weights = [0.5 0.5];
	end
	model.prior = [weights(1)*model1.prior(:); weights(2)*model2.prior(:)];
	model.mean = [model1.mean; model2.mean];
	model.icov = cat(1,model1.icov,model2.icov);
end


