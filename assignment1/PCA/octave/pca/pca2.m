function  [pdata,PC,V] = pca2(data)
% PCA2:   Perform principal component analysis using SVD.
%
% Usage:  [pdata,PC,V] = pca2(data)
%
% Principal Component Analysis (PCA) is a method of diagonlizing
% the covarince matrix of a data set. In real world applications,
% this means that the PCA can find the bases which maximize the
% variance and allow for easier dimesional reduction. There are
% several algorithms for implementing this analysis. This
% method is based on its important relationship with singular value
% decomposition.
%
% Arguments:
%     data - matrix (MxN) where M is the number of dimensions of
%            the data set and N is the number of data set trials.
%
%            In other words, each row is one type of variable and
%            each column is is observation of the data set (with M
%            variables)
%
% Outputs:
%  pdata - the projected data. This MxN matrix is the original data
%          set projected on to the principal component basis. This
%          is basically a rotation PC of the original basis.
%     PC - the principal components. Each column is a principal
%          component. The i'th principal component corresponds to
%          the i'th largest singular value V(i).
%      V - a vector of variances where V(i) corresponds to
%          PC(:,i). The variance quantifies the amount of 'spread'
%          of the data set along each dimension. It is solely a
%          second order statistic of the data.
%
%  24 August 2002
%  Jon Shlens | jonshlens [at] ucsd.edu
%


%%%%%%%%%%%%%%%%%%%%%%
% 1. FIND SIZE OF DATA
%%%%%%%%%%%%%%%%%%%%%%

[M,N] = size(data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. SUBTRACT MEAN OFF OF DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the mean for each dimension (variable)
mn =  mean(data,2);

% subtract off the mean for each dimension (variable)
data = data - repmat(mn,1,N);

%%%%%%%%%%%%%%%%%%%%
% 3. SVD DOES IT ALL
%%%%%%%%%%%%%%%%%%%%

[u,S,PC] = svd(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. CALCULATE THE VARIANCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Matlab function "svd" produces a diagonal matrix instead of a
% single vector. This can be fixed by a simple procedure.

S = diag(S);

% The singular values S from SVD are proportional to the square
% root of variances V. The proportionality constant is the
% 'unbiased' covariance estimator 1/(N-1).

V = 1 / (N-1) * S.*S;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. PROJECT THE DATA ON TO THE PC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just a rotation (or matrix multiplication) by the new basis

pdata = PC' * data;


