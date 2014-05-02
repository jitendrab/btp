% HMMCOMPILE
%
% Script to compile the mex file hmmlogsum.c. 
oldpath = pwd;
np = which('hmminit');
[newpath,fname] = fileparts(np);
cd(newpath);
mex hmmlogsum.c
cd(oldpath);
