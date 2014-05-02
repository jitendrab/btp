% Mixture of Gaussian Hidden Markov Model Toolbox
% Version 0.9.4 24-Jul-2012
%
%Hidden Markov Model functions
%-----------------------------
%hmminit            initialize HMM with MoG state models
%hmminitwithsegments  initialize HMM with time segments
%hmmp               HMM probabilities
%hmmlogp            HMM log-probabilities
%hmmem              EM algorithm to optimize HMM parameters
%hmmem_weighted     weighted EM algorithm to optimize HMM parameters
%hmmlogem           log-version of EM algorithm to optimize HMM
%hmmtimesegment     Segment a timeseries using an HMM model
%hmmviterbi         optimal path through a HMM model
%hmm2mdhmm          convert a standard to a minimum-duration HMM
%mogmodels2hmm      combine a set of MoG's to a HMM
%hmmcombinestates   combine states into a fused state
%hmmc               train a HMM per class
%hmmc2              train a HMM per class, where the train seq. are
%                   concatenated
%
%Mixture of Gaussian functions
%-----------------------------
%moginit            initialize MoG model
%mogp               estimate MoG probabilities
%moglogp            estimate MoG log-probabilities
%mogem              EM algorithm to optimize MoG parameters
%mog2prtools        convert to Prtools mapping
%mogcat             concatenate two MoG models
%mogfindbestpair    find the best MoG pair
%
%Dataset generation
%------------------
%gendatmoghmm       generate data from a HMM model
%gendatmog          generate data from a MoG model
%randmvn            generate multivariate normal (Gaussian) data
%
%Support functions
%-----------------
%hmmforwardbackward    standard forward-backward algorithm for HMM
%hmmlogforwardbackward the corresponding log-version
%mdhmmforwardbackward  forward-backward algorithm for HMM with minimum
%                   duration constraints
%mdhmmlogforwardbackward  the corresponding log-version
%unitnorm           normalize row/column of a matrix to sum=1
%vecunitnorm        normalize a vector to sum=1
%hmmlogsum          mex file for fast computation of a sum of logs
%hmmcompile         script to compile hmmlogsum
%
%Example functions
%-----------------
%tsthmm             generate data, train new HMM on data
%
% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

