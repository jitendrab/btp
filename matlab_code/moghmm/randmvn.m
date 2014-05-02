%RANDMVG Generate multivariate Normally distributed random numbers
%
%      X = RANDMVN(N,MN,SIGMA)
%
% Generates N objects from a MultiVariate Normal (Gaussian)
% Distribution, with a mean MN and a covariance matrix SIGMA.
%
% See also: randn

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function x = randmvn(N,mn,sigma)

%Parse input parameters
% Check the mean:
[h,dim] = size(mn);
if (dim==1)&&(h>1)
	dim = h;
	mn = mn';
end
% Check the cov-matrix
[h1,h2] = size(sigma);
if (h1~=h2)
	error('Covariance matrix should be square.');
end
if (h1~=dim)
	error('Sizes cov. matrix and mean vector do not match.');
end
if N < 1
    N=1;
end
N=fix(N); %Make sure N is a whole number

%Start of the program
x=randn(N,dim);  %generate the samples using built in Matlab function
xmean=mean(x); %calculate the sample mean

%subtract off the mean when N > 1
%This removes any mean from the Matlab generated numbers
%as N increases the Matlab mean approaches zero
if N>1
	x=x-repmat(xmean,[N,1]);
end

%Computes the Cholesky decomposition of the given variance
%It uses a different method depending on wheter or not variance is positive
%definite. This is used because the the variance is needed which is the
%square root of the covariance
[R,p]=chol(sigma);
if p>0
    x=x*sqrtm(sigma);
else
    x=x*R;
end

%Add back on the desired mean
x=x+repmat(mn,[N,1]);

    
