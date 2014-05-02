% generate sequences from two classes:
a = gendatseq([10 10],2,3);
mildisp(a)
scatterd(a)

% train a HMM classifier on this:
w = hmmc(a,1,[1 1],0.1,20);

% classify
a*w*labeld

