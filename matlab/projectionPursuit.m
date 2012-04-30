% [2] M = number of sources signals and signal mixtures
M = 2;
% [1e4] N = number of data points per signal.
N = 1e4;

%Load data, each M=2 rows contains a different source signal.  
% Each row has N columns (signal values.)

% Load standard matlab sounds (from MatLab's datafun directory)
% Set variance of each source to unity. 
load chirp; s1 = y(1:N) ; s1 = s1 - mean(s1); s1 = s1' / std(s1);
load gong; s2 = y(1:N); s2 = s2 - mean(s2); s2 = s2' / std(s2);

% Combine sources into vector variable s. 
s = [s1; s2];

% Make mixing matrix 
A = randn (M, M);

% Listen to source signals .... 
% [10000] Fs Sample rate of sound.
listen=true;
Fs = 10000;
if listen  soundsc(s(1,:), Fs); soundsc(s(2,:), Fs); end;

% Plot histogram of each source signal - 
% this approximates pdf of each source. 

figure(3); hist(s(1,:), 50); drawnow;
figure(4);hist(s(2,:), 50); drawnow;

% Make N mixtures x from M sources signals s. 
x = A *s;

% Sphere mixtures using SVD.
[U D V] = svd(x', 0);

% Set new x to be left singular vectors of old x.
z = U;

% Each eigenvector has unit length,
%	but we want unit variance mixtures .... 
z = z ./ repmat(std(z,1),N,1);

% Initialise unmixing vector to random vector;
w = rand(1, M)';
%  ... with unit length
w = w/norm(w);

% Initialize the first y
y = w' * z;

% Print out initial correlation between each estimated source y and every source signal s.
fprintf('Initial correlations of source extracted signals\n');
% rinitial= abs(r(M+1: 2 * M, 1:M));
r1=corrcoef([y;s1]');
r2=corrcoef([y;s2]');
rinitial=abs([r1(1,2) r2(1,2)])

maxiter =100;

% Make array hs to store values of function and gradient magnitude.
Ks=zero(maxiter, 1);
gs=zero(maxiter, 1);

%Begin gradient ascent on K ... 
% Define know optimal weight vector ... 
wopt=[-0.6125   0.794];

for iter=1:maxiter;
	% Get estimated source signal signal, y. 
	y = w' * z;

	% Get estimated kurtosis.
	K = mean (y .^4) - 3;

	% Find gradient  @K/@w ...
	y3 = y.^3;
	yy3 = repmat(y3,2,1);
	g = mean( (z.*yy3)')';

	% Update w to increase K ...
	w = w + eta *g;

	% Set length of w to unity ...
	w = w/ norm(w);

	% Record h and angle between wopt and gradient
	Ks(iter)= K;
	gs(iter) = subspace(g, wopt');

end;

% Plot change in K and gradient/ wopt angle during optimisation. 
jfig(1); plot(Ks, 'k');
title('Function values - Kurtosis');
xlabel('Iteration');
ylabel('K(y)');
jfig(2);plot(gs, 'k');
title('Angle \alpha Between g and Final Weight Vector');
xlabel('Iteration');ylabel('\alpha');

% Print out final correlations ... 
r = corrcoef([y; s]');
fprintf('Final correlations between source and extracted signals ...\n');
r1 = corrcoef([y;s1]');
r2 = corrcoef([y;s2]');
rfinal=abs([r1(1,2) r2(1,2)])

% Listen to extracted signal ...
if listen soundsc(y,Fs); end;