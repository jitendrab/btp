% Test the HMM code:

% Generate dataset:
dim = 2;
sd = 3;
hmmgen = hmminit([2 2],dim,sd);
nrsegm = 25;
minT = 10;
[x,labx] = gendatmoghmm(hmmgen,nrsegm,minT);

% standard settings:
M = [2 2];
dim = size(x,2);
rvar = 4;
stopcrit.maxiters = 50;
stopcrit.minllimpr = 1e-5;
reg = 0.01;

% Make a HMM and fit it using EM:
hmm0 = hmminit(M,dim,rvar);
%hmm0 = hmminitwithsegments(x,2);
%[p,p_k]=hmmp(x,hmm0);
%[g0,g10,sumxi0,ll0,g] = hmmforwardbackward(x,hmm0,p,p_k);
%[g1,g11,sumxi1,ll1,logg] = hmmlogforwardbackward(x,hmm0,p,p_k);
tic
hmm1 = hmmem(x,hmm0,stopcrit,reg);
toc
[pth,tmp,sortofprob] = hmmviterbi(x,hmm1);
%confmat(labx,pth)
tic
hmm2 = hmmlogem(x,hmm0,stopcrit,reg);
toc
% Find the most likeli path and show it:
pth = hmmviterbi(x,hmm2);
% show:
figure(1); clf; plot(x); hold on; plot(labx,'k');
%confmat(labx,pth), plot(pth,'r')

figure(2); clf; scatterd(x);
for i=1:length(hmm1.pdf)
	hold on; scatterd(hmm1.pdf{i}.mean,'ro');
end
figure(3); clf; scatterd(x);
for i=1:length(hmm1.pdf)
	hold on; scatterd(hmm2.pdf{i}.mean,'ro');
end
