function  [pdata,PC,V] = pca1(data)
% PCA1:   Perform principal component analysis using covariance.
%
% Usage:  [pdata,PC,V] = pca1(data)
%
% Principal Component Analysis (PCA) is a method of diagonlizing
% the covarince matrix of a data set. In real world applications,
% this means that the PCA can find the bases which maximize the
% variance and allow for easier dimesional reduction. There are
% several algorithms for implementing this analysis. This
% method involves the explicit calculation of the covariance.
%
% Arguments:
%
%     data - matrix (MxN) where M is the number of dimensions of
%            the data set and N is the number of data set trials.
%
%            In other words, each row is one type of variable and
%            each column is is observation of the data set (with M
%            variables)
%
% Outputs:
%
%  pdata - the projected data. This MxN matrix is the original data
%          set projected on to the principal component basis. This
%          is basically a rotation PC of the original basis.
%     PC - the principal components. Each column is a principal
%          component. The i'th principal component corresponds to
%          the i'th largest variance V(i).
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. CALCULATE THE COVARIANCE MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that algorithm is using 'unbiased' estimator by multplying
% by 1/(N-1) instead of 1/N.

covariance = 1 / (N-1) * data * data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. A COMPUTATIONAL COMMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is just a comment and no code.

% This algorithm is being explicit about all of the matrix
% computations necessary.  A more efficient method is to note that
% Matlab has a single function (commented below) that can supercede
% all of the above code. It both subtracts of the mean and
% calculates the covariance.

% covariance = cov(data');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. FIND THE EIGENVECTORS OF THE COVARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One can prove that a matrix is always diagonalized by a matrix of
% it eigenvectors. Thus, the eigenvectors diagonalize the
% covariance matrix and therefore are the principal components.

[PC, V] = eig(covariance);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. GRAB DIAGONAL OF THE EIGENVALUE MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Matllab function "eig" produces a diagonal matrix instead of
% a single vector. This can be fixed by a simple procedure.

V = diag(V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. SORT THE VARIANCES (EIGENVALUES) IN DESCENDING ORDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The Matlab function "sort" only sorts in ascending order. Because
% all of the singular values are positive. I can use a trick of
% computing sort(-1*V) in order to sort in descending order.

[junk, rindices] = sort(-1*V);

% Use the matrix index to sort V and PC
V  = V(rindices);
PC = PC(:,rindices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. PROJECT THE DATA INTO THE PC BASIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% project the original data set
pdata = PC' * data;
