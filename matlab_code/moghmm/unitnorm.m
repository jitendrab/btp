%UNITNORM Normalize a matrix
%
%     [A, SCALE] = UNITNORM(A,DIM)
%
% Rescale each row (DIM=2) or column (DIM=1) in matrix A, and return the
% new matrix with the scaling factor SCALE.
%
% See also: VECUNITNORM

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function [A, scale] = unitnorm(A,dim)

if nargin<2
	scale = sum(A(:));
	if scale==0
		scale=1;
	else
		A = A./scale;
	end

else
	if dim==2 %most occuring...
		scale = sum(A,2);
		scale(scale==0) = 1;
		A = A./repmat(scale,1,size(A,2));
	elseif dim==1
		scale = sum(A,1);
		scale(scale==0) = 1;
		A = A./repmat(scale,size(A,1),1);
	else
		error('dim should be 1 or 2.');
	end
end
