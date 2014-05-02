%UNITNORM Normalize a vector
%
%     [X, SCALE] = UNITNORM(X)
%
% Rescale the vector to unit norm, and return the SCALE.
%
% See also: unitnorm

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function [x, scale] = vecunitnorm(x)

scale = sum(x);
if scale==0
	scale = 1;
else
	x = x/scale;
end
