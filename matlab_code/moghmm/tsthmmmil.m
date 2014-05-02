% does the hmm mil work??
hmmgen = hmminit([1 2],2,4,0.1);

nrseq = 20;
nrsegm = 10;
MD = 1;
a = []; laba = []; baglab = [];
for i=1:nrseq
	[x,labx] = gendatmoghmm(hmmgen,20,6);
	a = [a;x];
	laba = [laba; labx];
	baglab = [baglab; repmat(i,size(x,1),1)];
end

x = genmil(a,genmillabels(laba,1),baglab);
w = hmmc(x,'presence',[2 2],0.01);
confmat(x*w)

scatterd(dataset(x,x*w*labeld));
