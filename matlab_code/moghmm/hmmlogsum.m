%HMMLOGSUM Computes the log of a sum of exponentials
%
%      out = hmmlogsum(a,dim)
%
%   out = log( sum_i( exp(a_i)) ) 

function out = hmmlogsum(a,dim)
if nargin<2
	dim = 1;
end

out = log(sum(exp(a),dim));

return
