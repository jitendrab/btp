%        A = GENDATSEQ(N,DIM,SD,NRSEGM)
%
% Generate N sequences using two (randomly generated) HMM's (N(1) from
% the positive HMM, and N(2) from the negative HMM). The HMMs generate
% data in a DIM dimensional space, with a spread of SD. Each HMM
% generates NRSEGM segments in each sequence.
%
function a = gendatseq(N,dim,sd,nrsegm)
if nargin<4
	nrsegm = 15;
end
if nargin<3 || isempty(sd)
	sd = 3;
end
if nargin<2 || isempty(dim)
	dim = 2;
end
if nargin<1 || isempty(N)
	N = [10 10];
end
if length(N)<2
	N = [N N];
end
% generate sequences from two classes:
minT = 1;
x = [];
baglab = [];
lab = [];

hmmgen_pos = hmminit([1 1],dim,sd);
hmmgen_neg = hmminit([1 1],dim,2*sd);
for i=1:N(1)
	xpos = gendatmoghmm(hmmgen_pos,nrsegm,minT);
	m = size(xpos,1);
	x = [x; xpos];
	baglab = [baglab; repmat(i,m,1)];
	lab = [lab; repmat('positive',m,1)];
end
for i=1:N(2)
	xneg = gendatmoghmm(hmmgen_neg,nrsegm,minT);
	m = size(xneg,1);
	x = [x; xneg];
	baglab = [baglab; repmat(i+N(1),m,1)];
	lab = [lab; repmat('negative',m,1)];
end

a = genmil(x,lab,baglab);

return
