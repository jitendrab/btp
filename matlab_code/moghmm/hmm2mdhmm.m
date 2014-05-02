%HMM2MDHMM Add minimum duration constraint to standard HMM
%
%        MODEL = HMM2MDHMM(MODEL,MD)
%
% Convert a standard HMM to a Minimum Duration HMM, with a minimum
% duration of MD.  Note that ALL states have the same minimal duration.
% If you don't like that, please program it yourself:-)

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function model = hmm2mdhmm(hmm,md)
if nargin<2
	md = [];
end

% Check:
if isempty(md) || (md<=1)
	warning('moghmm:minDurationOfOne',...
        'Minimum duration is 1, now using the standard HMM.');
	model = hmm;
	model.md = 1;
	return;
end

% initialize:
Q = length(hmm.prior);

% prior is easy:
model.prior = zeros(Q*md,1);
model.prior(1:md:Q*md,1) = hmm.prior;

% transition is harder:
% first create a 'diagonal' matrix where the diagonal is one element
% off:
z1 = zeros(Q*md,1);
z2 = z1; z2(2)=1;
trans = toeplitz(z1,z2);
% then copy the elements on the right spot:
for i=1:Q
	trans(i*md,i*md) = hmm.trans(i,i);
	for j=i+1:Q
		trans(i*md,(j-1)*md+1) = hmm.trans(i,j);
		trans(j*md,(i-1)*md+1) = hmm.trans(j,i);
	end
end
model.trans = sparse(trans);

% all the pdf's do not change;
model.pdf = hmm.pdf;

% finally add the fact that we are working with minimum duration:
model.md = md;

end


