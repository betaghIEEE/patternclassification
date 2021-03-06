applyW.m                                                                                            0100644 0017072 0006200 00000002231 07450203724 011144  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function [Out1, Out2] = applyW(mixedsig, W, ...)

# File:  applyW.m 
#
# Description: Script for applying the inverse mixing matrix W to a matrix of data
#
# Examples: 
# [icasig, meanvalues] = applyW(mixedsig, W)
# Gives the ica matrix and the mean values of mixedsig.  These means are not put 
# back into icasig.
#
# [icasig] = applyW(mixedsig, W, 'mean')
# Gives the ica matrix only, with the means put back into icasig.
#
# Created: February, 2002, Dan Ryan

if nargin==2
  [mixedsig, meanvalues] = remmean(mixedsig);
  icasig = W * mixedsig;
elseif (nargin<2)
  error('Not enough arguments');
  exit(-1);
elseif (nargin>3)
  error('Too many arguments');
  exit(-1);
elseif nargin==3
  param = va_arg();
  if strcmp('mean', param)
    icasig = W * mixedsig;
  else
    error(['Unrecognized parameter: ''' param '''']);
    exit(-1);
  endif
endif

# Determine what to output
if nargout == 1
  Out1 = icasig;
else
  Out1 = icasig;
  Out2 = meanvalues;
endif
# Plot to make sure they look reasonable:
icaplot('histogram', icasig);

# Say how big the ascii files are:
length = size(icasig, 2);
printf("The ascii files will have %i data points.\n", length);

endfunction





                                                                                                                                                                                                                                                                                                                                                                       demosig.m                                                                                           0100644 0017072 0006200 00000001315 07450203724 011321  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function [sig,mixedsig]=demosig();
%
% function [sig,mixedsig]=demosig();
% 
% Returns artificially generated test signals, sig, and mixed
% signals, mixedsig. Signals are row vectors of
% matrices. Input mixedsig to FastICA to see how it works.

% 2 Apr 1998 Aapo Hyv?rinen

%create source signals (independent components)
N=500; %data size

v=[0:N-1];
sig=[];
sig(1,:)=sin(v/2); %sinusoid
sig(2,:)=((rem(v,23)-11)/9).^5; %funny curve
sig(3,:)=((rem(v,27)-13)/9); %saw-tooth
sig(4,:)=((rand(1,N)<.5)*2-1).*log(rand(1,N)); %impulsive noise

for t=1:4
sig(t,:)=sig(t,:)/std(sig(t,:));
end

%remove mean (not really necessary)

[sig mean]=remmean(sig);

%create mixtures

Aorig=rand(size(sig,1));
mixedsig=(Aorig*sig);
                                                                                                                                                                                                                                                                                                                   dispsig.m                                                                                           0100644 0017072 0006200 00000000561 07450203724 011336  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function dispsig(signalMatrix, range, titlestr);
%DISPSIG - deprecated!
%
% Please use icaplot instead.
%
%   See also ICAPLOT

% 24.8.1998
% Hugo G?vert

fprintf('\nNote: DISPSIG is now deprecated! Please use ICAPLOT.\n');

if nargin < 3, titlestr = ''; end
if nargin < 2, range = 1:size(signalMatrix, 1); end

icaplot('dispsig',signalMatrix',0,range,range,titlestr);
                                                                                                                                               fastica.m                                                                                           0100644 0017072 0006200 00000044130 07450203724 011306  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function [Out1, Out2, Out3] = fastica(mixedsig, ...)
%FASTICA - Fast Independent Component Analysis
%
% FastICA for Matlab 5.x
% Version 2.1, January 15 2001
% Copyright (c) Hugo G?vert, Jarmo Hurri, Jaakko S?rel?, and Aapo Hyv?rinen.
%
% FASTICA(mixedsig) estimates the independent components from given
% multidimensional signals. Each row of matrix mixedsig is one
% observed signal.  FASTICA uses Hyvarinen's fixed-point algorithm,
% see http://www.cis.hut.fi/projects/ica/fastica/. Output from the
% function depends on the number output arguments:
%
% [icasig] = FASTICA (mixedsig); the rows of icasig contain the
% estimated independent components.
%
% [icasig, A, W] = FASTICA (mixedsig); outputs the estimated separating
% matrix W and the corresponding mixing matrix A.
%
% [A, W] = FASTICA (mixedsig); gives only the estimated mixing matrix
% A and the separating matrix W.
%
% Some optional arguments induce other output formats, see below.
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
% FASTICA can be called with numerous optional arguments. Optional
% arguments are given in parameter pairs, so that first argument is
% the name of the parameter and the next argument is the value for
% that parameter. Optional parameter pairs can be given in any order.
%
% OPTIONAL PARAMETERS:
%
% Parameter name        Values and description
%
%======================================================================
% --Basic parameters in fixed-point algorithm:
%
% 'approach'            (string) The decorrelation approach used. Can be
%                       symmetric ('symm'), i.e. estimate all the
%                       independent component in parallel, or
%                       deflation ('defl'), i.e. estimate independent
%                       component one-by-one like in projection pursuit.
%                       Default is 'defl'.
%
% 'numOfIC'             (integer) Number of independent components to
%                       be estimated. Default equals the dimension of data.
%
%======================================================================
% --Choosing the nonlinearity:
%
% 'g'                   (string) Chooses the nonlinearity g used in 
%                       the fixed-point algorithm. Possible values:
%
%                       Value of 'g':      Nonlinearity used:
%                       'pow3' (default)   g(u)=u^3
%                       'tanh'             g(u)=tanh(a1*u)
%                       'gauss             g(u)=u*exp(-a2*u^2/2)
%                       'skew'             g(u)=u^2
% 
% 'finetune'		(string) Chooses the nonlinearity g used when 
%                       fine-tuning. In addition to same values
%                       as for 'g', the possible value 'finetune' is:
%                       'off'              fine-tuning is disabled.
%
% 'a1'                  (number) Parameter a1 used when g='tanh'.
%                       Default is 1.
% 'a2'                  (number) Parameter a2 used when g='gaus'.
%                       Default is 1.
%
% 'mu'			(number) Step size. Default is 1.
%                       If the value of mu is other than 1, then the
%                       program will use the stabilized version of the
%                       algorithm (see also parameter 'stabilization').
%
%
% 'stabilization'       (string) Values 'on' or 'off'. Default 'off'. 
%                       This parameter controls wether the program uses
%                       the stabilized version of the algorithm or
%                       not. If the stabilization is on, then the value
%                       of mu can momentarily be halved if the program
%                       senses that the algorithm is stuck between two
%                       points (this is called a stroke). Also if there
%                       is no convergence before half of the maximum
%                       number of iterations has been reached then mu
%                       will be halved for the rest of the rounds.
% 
%======================================================================
% --Controlling convergence:
%
% 'epsilon'             (number) Stopping criterion. Default is 0.0001.
%
% 'maxNumIterations'    (integer) Maximum number of iterations.
%                       Default is 1000.
%
% 'maxFinetune'         (integer) Maximum number of iterations in 
%                       fine-tuning. Default 100.
%
% 'sampleSize'          (number) [0 - 1] Percentage of samples used in
%                       one iteration. Samples are chosen in random.
%                       Default is 1 (all samples).
%
% 'initGuess'           (matrix) Initial guess for A. Default is random.
%                       You can now do a "one more" like this: 
%                       [ica, A, W] = fastica(mix, 'numOfIC',3);
%                       [ica2, A2, W2] = fastica(mix, 'initGuess', A, 'numOfIC', 4);
%
%======================================================================
% --Graphics and text output:
%
% 'verbose'             (string) Either 'on' or 'off'. Default is
%                       'on': report progress of algorithm in text format.
%
% 'displayMode'         (string) Plot running estimates of independent
%                       components: 'signals', 'basis', 'filters' or
%                       'off'. Default is 'signals'.
%
% 'displayInterval'     Number of iterations between plots.
%                       Default is 1 (plot after every iteration).
%
%======================================================================
% --Controlling reduction of dimension and whitening:
%
% Reduction of dimension is controlled by 'firstEig' and 'lastEig', or
% alternatively by 'interactivePCA'. 
%
% 'firstEig'            (integer) This and 'lastEig' specify the range for
%                       eigenvalues that are retained, 'firstEig' is
%                       the index of largest eigenvalue to be
%                       retained. Default is 1.
%
% 'lastEig'             (integer) This is the index of the last (smallest)
%                       eigenvalue to be retained. Default equals the
%                       dimension of data.
%
% 'interactivePCA'      (string) Either 'on' or 'off'. When set 'on', the
%                       eigenvalues are shown to the user and the
%                       range can be specified interactively. Default
%                       is 'off'. Can also be set to 'gui'. Then the user
%                       can use the same GUI that's in FASTICAG.
%
% If you already know the eigenvalue decomposition of the covariance
% matrix, you can avoid computing it again by giving it with the
% following options:
%
% 'pcaE'                (matrix) Eigenvectors
% 'pcaD'                (matrix) Eigenvalues
%
% If you already know the whitened data, you can give it directly to
% the algorithm using the following options:
%
% 'whiteSig'            (matrix) Whitened signal
% 'whiteMat'            (matrix) Whitening matrix
% 'dewhiteMat'          (matrix) dewhitening matrix
%
% If values for all the 'whiteSig', 'whiteSig' and 'dewhiteMat' are
% suplied, they will be used in computing the ICA. PCA and whitening
% are not performed. Though 'mixedsig' is not used in the main
% algorithm it still must be entered - some values are still
% calculated from it.
%
% Performing preprocessing only is possible by the option:
%
% 'only'                (string) Compute only PCA i.e. reduction of
%                       dimension ('pca') or only PCA plus whitening
%                       ('white'). Default is 'all': do ICA estimation
%                       as well.  This option changes the output
%                       format accordingly. For example: 
%
%                       [whitesig, WM, DWM] = FASTICA(mixedsig, 
%                       'only', 'white') 
%                       returns the whitened signals, the whitening matrix
%                       (WM) and the dewhitening matrix (DWM). (See also
%                       WHITENV.) In FastICA the whitening matrix performs
%                       whitening and the reduction of dimension. Dewhitening
%                       matrix is the pseudoinverse of whitening matrix.
%                        
%                       [E, D] = FASTICA(mixedsig, 'only', 'pca') 
%                       returns the eigenvector (E) and diagonal 
%                       eigenvalue (D) matrices  containing the 
%                       selected subspaces. 
%
%======================================================================
% EXAMPLES
%
%       [icasig] = FASTICA (mixedsig, 'approach', 'symm', 'g', 'tanh');
%               Do ICA with tanh nonlinearity and in parallel (like
%               maximum likelihood estimation for supergaussian data).
%
%       [icasig] = FASTICA (mixedsig, 'lastEig', 10, 'numOfIC', 3);
%               Reduce dimension to 10, and estimate only 3
%               independent components.
%
%       [icasig] = FASTICA (mixedsig, 'verbose', 'off', 'displayMode', 'off');
%               Don't output convergence reports and don't plot
%               independent components.
%
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
%   See also FASTICAG

% 15.1.2001
% Hugo G?vert


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the mean and check the data

[mixedsig, mixedmean] = remmean(mixedsig);

[Dim, NumOfSampl] = size(mixedsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values for optional parameters

% All
verbose           = 'on';

% Default values for 'pcamat' parameters
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';

% Default values for 'fpica' parameters
approach          = 'defl';
numOfIC           = Dim;
g                 = 'pow3';
finetune          = 'off';
a1                = 1;
a2                = 1;
myy               = 1;
stabilization     = 'off';
epsilon           = 0.0001;
maxNumIterations  = 1000;
maxFinetune       = 5;
initState         = 'rand';
guess             = 0;
sampleSize        = 1;
displayMode       = 'signals';
displayInterval   = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for fastICA - i.e. this file

b_verbose = 1;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters passed to icaplot
% Added by DR, 12.18.2001

plottype = 'histogram';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the optional parameters

if(rem(nargin-1,2)==1)
     error('Optional parameters should always go by pairs');
     else
for i=1:(nargin-1)/2
% get the name and value of parameter
str_param = va_arg();
val_param = va_arg();

% change the value of parameter
if strcmp (str_param, 'stabilization')
     stabilization = val_param;
     elseif strcmp (str_param, 'maxFinetune')
     maxFinetune = val_param;
     elseif strcmp (str_param, 'sampleSize')
     sampleSize = val_param;
     elseif strcmp (str_param, 'verbose')
     verbose = val_param;
     % silence this program also
     if strcmp (verbose, 'off'), b_verbose = 0; end
     elseif strcmp (str_param, 'firstEig')
     firstEig = val_param;
     elseif strcmp (str_param, 'lastEig')
     lastEig = val_param;
     elseif strcmp (str_param, 'interactivePCA')
     interactivePCA = val_param;
     elseif strcmp (str_param, 'approach')
     approach = val_param;
     elseif strcmp (str_param, 'numOfIC')
     numOfIC = val_param;
     % User has suplied new value for numOfIC.
     % We'll use this information later on...
      userNumOfIC = 1;
      elseif strcmp (str_param, 'g')
      g = val_param;
      elseif strcmp (str_param, 'finetune')
      finetune = val_param;
      elseif strcmp (str_param, 'a1')
      a1 = val_param;
      elseif strcmp (str_param, 'a2')
      a2 = val_param;
      elseif strcmp (str_param, 'mu') | strcmp(str_param, 'myy')
      myy = val_param;
      elseif strcmp (str_param, 'epsilon')
      epsilon = val_param;
      elseif strcmp (str_param, 'maxNumIterations')
      maxNumIterations = val_param;
      elseif strcmp (str_param, 'initGuess')
      % no use setting 'guess' if the 'initState' is not set
      initState = 'guess';
      guess = val_param;
      elseif strcmp (str_param, 'displayMode')
      displayMode = val_param;
      elseif strcmp (str_param, 'displayInterval')
      displayInterval = val_param;
      elseif strcmp (str_param, 'pcaE')
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      E = val_param;
      elseif strcmp (str_param, 'pcaD')
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      D = val_param;
      elseif strcmp (str_param, 'whiteSig')
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whitesig = val_param;
      elseif strcmp (str_param, 'whiteMat')
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whiteningMatrix = val_param;
      elseif strcmp (str_param, 'dewhiteMat')
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      dewhiteningMatrix = val_param;
     elseif strcmp (str_param,  'only')
      % if the user only wants to calculate PCA or...
        if strcmp (val_param, 'pca')
	only = 1;
       elseif strcmp(val_param, 'white')
	only = 2;
       elseif strcmp(val_param, 'all')
	only = 3;
      end
% Added by DR, define the type of plot 
     elseif strcmp (str_param, 'plottype')
     plottype = val_param;
      
     else
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' str_param '''']);
    end;
  end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print information about data
if b_verbose
  fprintf('Number of signals: %d\n', Dim);
  fprintf('Number of samples: %d\n', NumOfSampl);
end

% Check if the data has been entered the wrong way,
% but warn only... it may be on purpose

if Dim > NumOfSampl
  if b_verbose
    fprintf('Warning: ');
    fprintf('The signal matrix may be oriented in the wrong way.\n');
    fprintf('In that case transpose the matrix.\n\n');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating PCA

% We need the results of PCA for whitening, but if we don't
% need to do whitening... then we dont need PCA...
if jumpWhitening == 3
  if b_verbose,
    fprintf ('Whitened signal and corresponding matrices supplied.\n');
    fprintf ('PCA calculations not needed.\n');
  end;
else
  
  % OK, so first we need to calculate PCA
  % Check to see if we already have the PCA data
  if jumpPCA == 2,
    if b_verbose,
      fprintf ('Values for PCA calculations suplied.\n');
      fprintf ('PCA calculations not needed.\n');
    end;
  else
    % display notice if the user entered one, but not both, of E and D.
    if (jumpPCA > 0) & (b_verbose),
      fprintf ('You must supply all of these in order to jump PCA:\n');
      fprintf ('''pcaE'', ''pcaD''.\n');
    end;
    
    % Calculate PCA
    [E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);
  end
end

% skip the rest if user only wanted PCA
if only > 1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Whitening the data
  
  % Check to see if the whitening is needed...
  if jumpWhitening == 3,
    if b_verbose,
      fprintf ('Whitening not needed.\n');
    end;
  else
    
    % Whitening is needed
    % display notice if the user entered some of the whitening info, but not all.
    if (jumpWhitening > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump whitening:\n');
      fprintf ('''whiteSig'', ''whiteMat'', ''dewhiteMat''.\n');
    end;
    
    % Calculate the whitening
    [whitesig, whiteningMatrix, dewhiteningMatrix] = whitenv ...
						     (mixedsig, E, D, verbose);
  end
  
end % if only > 1

% skip the rest if user only wanted PCA and whitening
if only > 2
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculating the ICA
  
  % Check some parameters
  % The dimension of the data may have been reduced during PCA calculations.
  % The original dimension is calculated from the data by default, and the
  % number of IC is by default set to equal that dimension.
  
  Dim = size(whitesig, 1);
  
  % The number of IC's must be less or equal to the dimension of data
  if numOfIC > Dim
    numOfIC = Dim;
    % Show warning only if verbose = 'on' and user suplied a value for 'numOfIC'
    if (b_verbose & userNumOfIC)
      fprintf('Warning: estimating only %d independent components\n', numOfIC);
      fprintf('(Can''t estimate more independent components than dimension of data)\n');
    end
  end
  
  % Calculate the ICA with fixed point algorithm.
  [A, W] = fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, approach, ...
                  numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon, ...
                  maxNumIterations, maxFinetune, initState, guess, sampleSize, ...
                  displayMode, displayInterval, verbose, plottype);
  
  % Check for valid return
  if ~isempty(W)
    % Add the mean back in.
    if b_verbose
      fprintf('Adding the mean back to the data.\n');
    end
    %icasig = W * mixedsig + (W * mixedmean) * ones(1, NumOfSampl);
    icasig = W * mixedsig;
    if b_verbose & ...
	  (max(abs(W * mixedmean)) > 1e-9) & ...
	  (strcmp(displayMode,'signals') | strcmp(displayMode,'on'))
      fprintf('Note that the plots don''t have the mean added.\n');
    end
  else
    icasig = [];
  end

end % if only > 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output depends on the number of output parameters
% and the 'only' parameter.

if only == 1    % only PCA
  Out1 = E;
  Out2 = D;
elseif only == 2  % only PCA & whitening
  if nargout == 2
    Out1 = whiteningMatrix;
    Out2 = dewhiteningMatrix;
  else
    Out1 = whitesig;
    Out2 = whiteningMatrix;
    Out3 = dewhiteningMatrix;
  end
else      % ICA
  if nargout == 2
    Out1 = A;
    Out2 = W;
  else
    Out1 = icasig;
    Out2 = A;
    Out3 = W;
  end
end

endfunction
                                                                                                                                                                                                                                                                                                                                                                                                                                        fpica.m                                                                                             0100644 0017072 0006200 00000063360 07450203724 010764  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function [A, W] = fpica(X, whiteningMatrix, dewhiteningMatrix, approach, ...
			numOfIC, g, finetune, a1, a2, myy, stabilization, ...
			epsilon, maxNumIterations, maxFinetune, initState, ...
			guess, sampleSize, displayMode, displayInterval, ...
			s_verbose, plottype);

%FPICA - Fixed point ICA. Main algorithm of FASTICA.
%
% [A, W] = fpica(whitesig, whiteningMatrix, dewhiteningMatrix, approach,
%        numOfIC, g, finetune, a1, a2, mu, stabilization, epsilon, 
%        maxNumIterations, maxFinetune, initState, guess, sampleSize,
%        displayMode, displayInterval, verbose);
% 
% Perform independent component analysis using Hyvarinen's fixed point
% algorithm. Outputs an estimate of the mixing matrix A and its inverse W.
%
% whitesig                              :the whitened data as row vectors
% whiteningMatrix                       :can be obtained with function whitenv
% dewhiteningMatrix                     :can be obtained with function whitenv
% approach      [ 'symm' | 'defl' ]     :the approach used (deflation or symmetric)
% numOfIC       [ 0 - Dim of whitesig ] :number of independent components estimated
% g             [ 'pow3' | 'tanh' |     :the nonlinearity used
%                 'gaus' | 'skew' ]     
% finetune      [same as g + 'off']     :the nonlinearity used in finetuning.
% a1                                    :parameter for tuning 'tanh'
% a2                                    :parameter for tuning 'gaus'
% mu                                    :step size in stabilized algorithm
% stabilization [ 'on' | 'off' ]        :if mu < 1 then automatically on
% epsilon                               :stopping criterion
% maxNumIterations                      :maximum number of iterations 
% maxFinetune                           :maximum number of iteretions for finetuning
% initState     [ 'rand' | 'guess' ]    :initial guess or random initial state. See below
% guess                                 :initial guess for A. Ignored if initState = 'rand'
% sampleSize    [ 0 - 1 ]               :percentage of the samples used in one iteration
% displayMode   [ 'signals' | 'basis' | :plot running estimate
%                 'filters' | 'off' ]
% displayInterval                       :number of iterations we take between plots
% verbose       [ 'on' | 'off' ]        :report progress in text format
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%       [nv, wm, dwm] = whitenv(vectors, E, D);
%       [A, W] = fpica(nv, wm, dwm);
%
%
% This function is needed by FASTICA and FASTICAG
%
%   See also FASTICA, FASTICAG, WHITENV, PCAMAT

% 15.1.2001
% Hugo G?vert

% Added by DR
% plottype  [ 'dispsig' | 'classic' |   ''  | 'complot'
%             'scatter' | 'compare' | 'sum' | 'sumerror' ] 
%
%                                       :mode of the output plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variable for stopping the ICA calculations from the GUI
global g_FastICA_interrupt=[];
if isempty(g_FastICA_interrupt)
  clear global g_FastICA_interrupt;
  interruptible = 0;
else
  interruptible = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values

if nargin < 3, error('Not enough arguments!'); end
[vectorSize, numSamples] = size(X);
if nargin < 20, s_verbose = 'on'; end
if nargin < 19, displayInterval = 1; end
if nargin < 18, displayMode = 'on'; end
if nargin < 17, sampleSize = 1; end
if nargin < 16, guess = 1; end
if nargin < 15, initState = 'rand'; end
if nargin < 14, maxFinetune = 100; end
if nargin < 13, maxNumIterations = 1000; end
if nargin < 12, epsilon = 0.0001; end
if nargin < 11, stabilization = 'on'; end
if nargin < 10, myy = 1; end
if nargin < 9, a2 = 1; end
if nargin < 8, a1 = 1; end
if nargin < 7, finetune = 'off'; end
if nargin < 6, g = 'pow3'; end
if nargin < 5, numOfIC = vectorSize; end     % vectorSize = Dim
if nargin < 4, approach = 'defl'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the data

if imag(X) <> 0
  error('Input has an imaginary part.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for verbose

str_param = lower(s_verbose);
     
if strcmp (str_param, 'on'),
   b_verbose = 1;
elseif strcmp (str_param,  'off'),
  b_verbose = 0;
else
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for approach

str_param = lower(approach);
if strcmp (str_param, 'symm'),
  approachMode = 1;
elseif strcmp (str_param, 'defl'),
  approachMode = 2;
else
  error(sprintf('Illegal value [ %s ] for parameter: ''approach''\n', approach));
end
if b_verbose, fprintf('Used approach [ %s ].\n', approach); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for numOfIC

if vectorSize < numOfIC
  error('Must have numOfIC <= Dimension!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the sampleSize
if sampleSize > 1
  sampleSize = 1;
  if b_verbose
    fprintf('Warning: Setting ''sampleSize'' to 1.\n');
  end  
elseif sampleSize < 1
  if (sampleSize * numSamples) < 1000
    sampleSize = min(1000/numSamples, 1);
    if b_verbose
      fprintf('Warning: Setting ''sampleSize'' to %0.3f (%d samples).\n', ...
	      sampleSize, floor(sampleSize * numSamples));
    end  
  end
end
if b_verbose
  if  b_verbose & (sampleSize < 1)
    fprintf('Using about %0.0f%% of the samples in random order in every step.\n',sampleSize*100);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for nonlinearity.

str_param = lower(g);
if strcmp (str_param, 'pow3'),
  gOrig = 10;
elseif strcmp (str_param, 'tanh'),
  gOrig = 20;
elseif strcmp (str_param, 'gaus') | strcmp (str_param, 'gauss'),
  gOrig = 30;
elseif strcmp (str_param, 'skew'),
  gOrig = 40;
else
  error(sprintf('Illegal value [ %s ] for parameter: ''g''\n', g));
endif

if sampleSize ~= 1
  gOrig = gOrig + 2;
end
if myy ~= 1
  gOrig = gOrig + 1;
end

if b_verbose,
  fprintf('Used nonlinearity [ %s ].\n', g);
end

finetuningEnabled = 1;
str_param = lower(finetune);
if strcmp (str_param, 'pow3'),
    gFine = 10 + 1;
elseif strcmp(str_param, 'tanh'),
     gFine = 20 + 1;
elseif strcmp(str_param, 'gaus'), 
    gFine = 30 + 1;
elseif strcmp(str_param, 'gauss'),
    gFine = 30 + 1;
elseif strcmp(str_param, 'skew'),
    gFine = 40 + 1;
elseif strcmp(str_param, 'off'),
    if (myy ~= 1),
     gFine = gOrig;
     else 
     gFine = gOrig + 1;
    endif
finetuningEnabled = 0;
else
error(sprintf('Illegal value [ %s ] for parameter: ''finetune''\n', ...
    finetune));
endif

if b_verbose & finetuningEnabled
  fprintf('Finetuning enabled (nonlinearity: [ %s ]).\n', finetune);
endif

str_param = lower(stabilization);
if strcmp (str_param, 'on'),
    stabilizationEnabled = 1;
elseif strcmp (str_param, 'off'),
    if myy ~= 1,
    stabilizationEnabled = 1;
    else
    stabilizationEnabled = 0;
    endif
else
  error(sprintf('Illegal value [ %s ] for parameter: ''stabilization''\n', ...
		stabilization)); 
endif

if b_verbose & stabilizationEnabled
  fprintf('Using stabilized algorithm.\n');
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some other parameters
myyOrig = myy;
% When we start fine-tuning we'll set myy = myyK * myy
myyK = 0.01;
% How many times do we try for convergence until we give up.
failureLimit = 5;


usedNlinearity = gOrig;
stroke = 0;
notFine = 1;
long = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for initial state.

str_param = lower(initState);
if strcmp (str_param, 'rand'),
  initialStateMode = 0;
elseif strcmp (str_param,'guess'),
    if size(guess,1) ~= size(whiteningMatrix,2),
     initialStateMode = 0;
     if b_verbose
      fprintf('Warning: size of initial guess is incorrect. Using random initial guess.\n');
     endif
    else
    initialStateMode = 1;
    if size(guess,2) < numOfIC
      if b_verbose
	fprintf('Warning: initial guess only for first %d components. Using random initial guess for others.\n', size(guess,2)); 
      end
      guess(:, size(guess, 2) + 1:numOfIC) = ...
					     rand(vectorSize,numOfIC-size(guess,2))-.5;
    elseif size(guess,2)>numOfIC
      guess=guess(:,1:numOfIC);
      fprintf('Warning: Initial guess too large. The excess column are dropped.\n');
    end
    if b_verbose, fprintf('Using initial guess.\n'); end
  end
else
  error(sprintf('Illegal value [ %s ] for parameter: ''initState''\n', initState));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for display mode.

str_param = lower(displayMode);
if strcmp (str_param, 'off') | strcmp (str_param, 'none'),
  usedDisplay = 0;
elseif strcmp (str_param, 'on') | strcmp (str_param, 'signals');
  usedDisplay = 1;
  if (b_verbose & (numSamples > 10000))
    fprintf('Warning: Data vectors are very long. Plotting may take long time.\n');
  end
  if (b_verbose & (numOfIC > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
elseif strcmp (str_param, 'basis'),
  usedDisplay = 2;
  if (b_verbose & (numOfIC > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
elseif strcmp (str_param, 'filters'),
  usedDisplay = 3;
  if (b_verbose & (vectorSize > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
else
  error(sprintf('Illegal value [ %s ] for parameter: ''displayMode''\n', displayMode));
end

% The displayInterval can't be less than 1...
if displayInterval < 1
  displayInterval = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b_verbose, fprintf('Starting ICA calculation...\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYMMETRIC APPROACH
if approachMode == 1,

  % set some parameters more...
  usedNlinearity = gOrig;
  stroke = 0;
  notFine = 1;
  long = 0;
  
  A = zeros(vectorSize, numOfIC);  % Dewhitened basis vectors.
  if initialStateMode == 0
    % Take random orthonormal initial vectors.
    B = orth(rand(vectorSize, numOfIC) - .5);
  elseif initialStateMode == 1
    % Use the given initial vector as the initial state
    B = whiteningMatrix * guess;
  end
  
  BOld = zeros(size(B));
  BOld2 = zeros(size(B));
  
  % This is the actual fixed-point iteration loop.
  for round = 1:maxNumIterations + 1,
    if round == maxNumIterations + 1,
      if b_verbose 
        fprintf('No convergence after %d steps\n', maxNumIterations);
        fprintf('Note that the plots are probably wrong.\n');
      end
      A=[];
      W=[];
      return;
    end
    
    if (interruptible & g_FastICA_interrupt)
      if b_verbose 
        fprintf('\n\nCalculation interrupted by the user\n');
      end
      if ~isempty(B)
	W = B' * whiteningMatrix;
	A = dewhiteningMatrix * B;
      else
	W = [];
	A = [];
      end
      return;
    end
    
    
    % Symmetric orthogonalization.
    B = B * real(inv(B' * B)^(1/2));
    
    % Test for termination condition. Note that we consider opposite
    % directions here as well.
    minAbsCos = min(abs(diag(B' * BOld)));
    minAbsCos2 = min(abs(diag(B' * BOld2)));
    
    if (1 - minAbsCos < epsilon)
      if finetuningEnabled & notFine
        if b_verbose, fprintf('Initial convergence, fine-tuning: \n'); end;
        notFine = 0;
        usedNlinearity = gFine;
        myy = myyK * myyOrig;
        BOld = zeros(size(B));
        BOld2 = zeros(size(B));
	
      else
        if b_verbose, fprintf('Convergence after %d steps\n', round); end
	
        % Calculate the de-whitened vectors.
        A = dewhiteningMatrix * B;
        break;
      end
    elseif stabilizationEnabled
      if (~stroke) & (1 - minAbsCos2 < epsilon)
	if b_verbose, fprintf('Stroke!\n'); end;
	stroke = myy;
	myy = .5*myy;
% "mod" doesn't exist in Octave, use this function explicitely:
	if usedNlinearity-2.*floor(usedNlinearity./2) == 0
	  usedNlinearity = usedNlinearity + 1;
	end
      elseif stroke
	myy = stroke;
	stroke = 0;
% "mod" doesn't exist in Octave, use this function explicitely:
	if (myy == 1) & (usedNlinearity-2.*floor(usedNlinearity./2) ~= 0)
	  usedNlinearity = usedNlinearity - 1;
	end
      elseif (~long) & (round>maxNumIterations/2)
	if b_verbose, fprintf('Taking long (reducing step size)\n'); end;
	long = 1;
	myy = .5*myy;
	if usedNlinearity-2.*floor(usedNlinearity./2) == 0
	  usedNlinearity = usedNlinearity + 1;
	end
      end
    end
    
    BOld2 = BOld;
    BOld = B;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show the progress...
    if b_verbose
      if round == 1
        fprintf('Step no. %d\n', round);
      else
        fprintf('Step no. %d, change in value of estimate: %.3f \n', round, 1 - minAbsCos);
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Also plot the current state...
    switch usedDisplay
     case 1
      if rem(round, displayInterval) == 0,
	% There was and may still be other displaymodes...
	% 1D signals
	icaplot(plottype,(X'*B)');
      end
     case 2
      if rem(round, displayInterval) == 0,
	% ... and now there are :-)
	% 1D basis
	A = dewhiteningMatrix * B;
	icaplot(plottype,A');
      end
     case 3
      if rem(round, displayInterval) == 0,
	% ... and now there are :-)
	% 1D filters
	W = B' * whiteningMatrix;
	icaplot(plottype,W);
      end
     otherwise
    end
    
    %drawnow;
    
    switch usedNlinearity
      % pow3
     case 10
      B = (X * (( X' * B) .^ 3)) / numSamples - 3 * B;
     case 11
      % optimoitu - epsilonin kokoisia eroja
      % t?m? on optimoitu koodi, katso vanha koodi esim.
      % aikaisemmista versioista kuten 2.0 beta3
      Y = X' * B;
      Gpow3 = Y .^ 3;
      Beta = sum(Y .* Gpow3);
      D = diag(1 ./ (Beta - 3 * numSamples));
      B = B + myy * B * (Y' * Gpow3 - diag(Beta)) * D;
     case 12
      Xsub=X(:, getSamples(numSamples, sampleSize));
      B = (Xsub * (( Xsub' * B) .^ 3)) / size(Xsub,2) - 3 * B;
     case 13
      % Optimoitu
      Ysub=X(:, getSamples(numSamples, sampleSize))' * B;
      Gpow3 = Ysub .^ 3;
      Beta = sum(Ysub .* Gpow3);
      D = diag(1 ./ (Beta - 3 * size(Ysub', 2)));
      B = B + myy * B * (Ysub' * Gpow3 - diag(Beta)) * D;
      
      % tanh
     case 20
      hypTan = tanh(a1 * X' * B);
      B = X * hypTan / numSamples - ...
	  ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / numSamples * ...
	  a1;
     case 21
      % optimoitu - epsilonin kokoisia 
      Y = X' * B;
      hypTan = tanh(a1 * Y);
      Beta = sum(Y .* hypTan);
      D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2)));
      B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;
     case 22
      Xsub=X(:, getSamples(numSamples, sampleSize));
      hypTan = tanh(a1 * Xsub' * B);
      B = Xsub * hypTan / size(Xsub, 2) - ...
	  ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / size(Xsub, 2) * a1;
     case 23
      % Optimoitu
      Y = X(:, getSamples(numSamples, sampleSize))' * B;
      hypTan = tanh(a1 * Y);
      Beta = sum(Y .* hypTan);
      D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2)));
      B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;
      
      % gauss
     case 30
      U = X' * B;
      Usquared=U .^ 2;
      ex = exp(-a2 * Usquared / 2);
      gauss =  U .* ex;
      dGauss = (1 - a2 * Usquared) .*ex;
      B = X * gauss / numSamples - ...
	  ones(size(B,1),1) * sum(dGauss)...
	  .* B / numSamples ;
     case 31
      % optimoitu
      Y = X' * B;
      ex = exp(-a2 * (Y .^ 2) / 2);
      gauss = Y .* ex;
      Beta = sum(Y .* gauss);
      D = diag(1 ./ (Beta - sum((1 - a2 * (Y .^ 2)) .* ex)));
      B = B + myy * B * (Y' * gauss - diag(Beta)) * D;
     case 32
      Xsub=X(:, getSamples(numSamples, sampleSize));
      U = Xsub' * B;
      Usquared=U .^ 2;
      ex = exp(-a2 * Usquared / 2);
      gauss =  U .* ex;
      dGauss = (1 - a2 * Usquared) .*ex;
      B = Xsub * gauss / size(Xsub,2) - ...
	  ones(size(B,1),1) * sum(dGauss)...
	  .* B / size(Xsub,2) ;
     case 33
      % Optimoitu
      Y = X(:, getSamples(numSamples, sampleSize))' * B;
      ex = exp(-a2 * (Y .^ 2) / 2);
      gauss = Y .* ex;
      Beta = sum(Y .* gauss);
      D = diag(1 ./ (Beta - sum((1 - a2 * (Y .^ 2)) .* ex)));
      B = B + myy * B * (Y' * gauss - diag(Beta)) * D;
      
      % skew
     case 40
      B = (X * ((X' * B) .^ 2)) / numSamples;
     case 41
      % Optimoitu
      Y = X' * B;
      Gskew = Y .^ 2;
      Beta = sum(Y .* Gskew);
      D = diag(1 ./ (Beta));
      B = B + myy * B * (Y' * Gskew - diag(Beta)) * D;
     case 42
      Xsub=X(:, getSamples(numSamples, sampleSize));
      B = (Xsub * ((Xsub' * B) .^ 2)) / size(Xsub,2);
     case 43
      % Uusi optimoitu
      Y = X(:, getSamples(numSamples, sampleSize))' * B;
      Gskew = Y .^ 2;
      Beta = sum(Y .* Gskew);
      D = diag(1 ./ (Beta));
      B = B + myy * B * (Y' * Gskew - diag(Beta)) * D;

     otherwise
      error('Code for desired nonlinearity not found!');
    end
  end

  
  % Calculate ICA filters.
  W = B' * whiteningMatrix;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Also plot the last one...
  switch usedDisplay
   case 1 
    % There was and may still be other displaymodes...
    % 1D signals
    icaplot(plottype,(X'*B)');
    %drawnow;
   case 2
    % ... and now there are :-)
    % 1D basis
    icaplot(plottype,A');
    %drawnow;
   case 3
    % ... and now there are :-)
    % 1D filters
    icaplot(plottype,W);
    %drawnow;
   otherwise
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFLATION APPROACH
if approachMode == 2
  
  B = zeros(vectorSize);
  
  % The search for a basis vector is repeated numOfIC times.
  round = 1;
  
  numFailures = 0;
  
  while round <= numOfIC,
    myy = myyOrig;
    usedNlinearity = gOrig;
    stroke = 0;
    notFine = 1;
    long = 0;
    endFinetuning = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show the progress...
    if b_verbose, fprintf('IC %d ', round); end
    
    % Take a random initial vector of lenght 1 and orthogonalize it
    % with respect to the other vectors.
    if initialStateMode == 0
      w = rand(vectorSize, 1) - .5;
    elseif initialStateMode == 1
      w=whiteningMatrix*guess(:,round);
    end
    w = w - B * B' * w;
    w = w / norm(w);
    
    wOld = zeros(size(w));
    wOld2 = zeros(size(w));
    
    % This is the actual fixed-point iteration loop.
    %    for i = 1 : maxNumIterations + 1
    i = 1;
    gabba = 1;
    while i <= maxNumIterations + gabba
      %drawnow;
% DR also removed from next line " & g_FastICA_interrupt
      if (interruptible)
        if b_verbose 
          fprintf('\n\nCalculation interrupted by the user\n');
        end
        return;
      end
      
      % Project the vector into the space orthogonal to the space
      % spanned by the earlier found basis vectors. Note that we can do
      % the projection with matrix B, since the zero entries do not
      % contribute to the projection.
      w = w - B * B' * w;
      w = w / norm(w);
      
      if notFine
	if i == maxNumIterations + 1
	  if b_verbose
	    fprintf('\nComponent number %d did not converge in %d iterations.\n', round, maxNumIterations);
	  end
	  round = round - 1;
	  numFailures = numFailures + 1;
	  if numFailures > failureLimit
	    if b_verbose
	      fprintf('Too many failures to converge (%d). Giving up.\n', numFailures);
	    end
	    if round == 0
	      A=[];
	      W=[];
	    end
	    return;
	  end
	  % numFailures > failurelimit
	  break;
	end
	% i == maxNumIterations + 1
      else
	% if notFine
	if i >= endFinetuning
	  wOld = w; % So the algorithm will stop on the next test...
	end
      end
      % if notFine
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Show the progress...
      if b_verbose, fprintf('.'); end;
      
      
      % Test for termination condition. Note that the algorithm has
      % converged if the direction of w and wOld is the same, this
      % is why we test the two cases.
      if norm(w - wOld) < epsilon | norm(w + wOld) < epsilon
        if finetuningEnabled & notFine
          if b_verbose, fprintf('Initial convergence, fine-tuning: '); end;
          notFine = 0;
	  gabba = maxFinetune;
          wOld = zeros(size(w));
          wOld2 = zeros(size(w));
          usedNlinearity = gFine;
          myy = myyK * myyOrig;
	  
	  endFinetuning = maxFinetune + i;
	  
        else
          numFailures = 0;
          % Save the vector
          B(:, round) = w;
	  
          % Calculate the de-whitened vector.
          A(:,round) = dewhiteningMatrix * w;
          % Calculate ICA filter.
          W(round,:) = w' * whiteningMatrix;
	  
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Show the progress...
          if b_verbose, fprintf('computed ( %d steps ) \n', i); end
	  
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Also plot the current state...
          switch usedDisplay
	   case 1
	    if rem(round, displayInterval) == 0,
	      % There was and may still be other displaymodes...   
	      % 1D signals
	      temp = X'*B;
	      icaplot(plottype,temp(:,1:numOfIC)');
	      %drawnow;
	    end
	   case 2
	    if rem(round, displayInterval) == 0,
	      % ... and now there are :-) 
	      % 1D basis
	      icaplot(plottype,A');
	      %drawnow;
	    end
	   case 3
	    if rem(round, displayInterval) == 0,
	      % ... and now there are :-) 
	      % 1D filters
	      icaplot(plottype,W);
	      %drawnow;
	    end
          end
	  % switch usedDisplay
	  break; % IC ready - next...
        end
	%if finetuningEnabled & notFine
      elseif stabilizationEnabled
	if (~stroke) & (norm(w - wOld2) < epsilon | norm(w + wOld2) < ...
			epsilon)
	  stroke = myy;
	  if b_verbose, fprintf('Stroke!'); end;
	  myy = .5*myy;
	  if mod(usedNlinearity,2) == 0
	    usedNlinearity = usedNlinearity + 1;
	  end
	elseif stroke
	  myy = stroke;
	  stroke = 0;
	  if (myy == 1) & (mod(usedNlinearity,2) ~= 0)
	    usedNlinearity = usedNlinearity - 1;
	  end
	elseif (notFine) & (~long) & (i > maxNumIterations / 2)
	  if b_verbose, fprintf('Taking long (reducing step size) '); end;
	  long = 1;
	  myy = .5*myy;
	  if mod(usedNlinearity,2) == 0
	    usedNlinearity = usedNlinearity + 1;
	  end
	end
      end
      
      wOld2 = wOld;
      wOld = w;
      
      switch usedNlinearity
	% pow3
       case 10
	w = (X * ((X' * w) .^ 3)) / numSamples - 3 * w;
       case 11
	EXGpow3 = (X * ((X' * w) .^ 3)) / numSamples;
	Beta = w' * EXGpow3;
	w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
       case 12
	Xsub=X(:,getSamples(numSamples, sampleSize));
	w = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2) - 3 * w;
       case 13
	Xsub=X(:,getSamples(numSamples, sampleSize));
	EXGpow3 = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2);
	Beta = w' * EXGpow3;
	w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
	% tanh
       case 20
	hypTan = tanh(a1 * X' * w);
	w = (X * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / numSamples;
       case 21
	hypTan = tanh(a1 * X' * w);
	Beta = w' * X * hypTan;
	w = w - myy * ((X * hypTan - Beta * w) / ...
		       (a1 * sum((1-hypTan .^2)') - Beta));
       case 22
	Xsub=X(:,getSamples(numSamples, sampleSize));
	hypTan = tanh(a1 * Xsub' * w);
	w = (Xsub * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / size(Xsub, 2);
       case 23
	Xsub=X(:,getSamples(numSamples, sampleSize));
	hypTan = tanh(a1 * Xsub' * w);
	Beta = w' * Xsub * hypTan;
	w = w - myy * ((Xsub * hypTan - Beta * w) / ...
		       (a1 * sum((1-hypTan .^2)') - Beta));
	% gauss
       case 30
	% This has been split for performance reasons.
	u = X' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	w = (X * gauss - sum(dGauss)' * w) / numSamples;
       case 31
	u = X' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	Beta = w' * X * gauss;
	w = w - myy * ((X * gauss - Beta * w) / ...
		       (sum(dGauss)' - Beta));
       case 32
	Xsub=X(:,getSamples(numSamples, sampleSize));
	u = Xsub' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	w = (Xsub * gauss - sum(dGauss)' * w) / size(Xsub, 2);
       case 33
	Xsub=X(:,getSamples(numSamples, sampleSize));
	u = Xsub' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	Beta = w' * Xsub * gauss;
	w = w - myy * ((Xsub * gauss - Beta * w) / ...
		       (sum(dGauss)' - Beta));
	% skew
       case 40
	w = (X * ((X' * w) .^ 2)) / numSamples;
       case 41
	EXGskew = (X * ((X' * w) .^ 2)) / numSamples;
	Beta = w' * EXGskew;
	w = w - myy * (EXGskew - Beta*w)/(-Beta);
       case 42
	Xsub=X(:,getSamples(numSamples, sampleSize));
	w = (Xub * ((Xub' * w) .^ 2)) / size(Xsub, 2);
       case 43
	Xsub=X(:,getSamples(numSamples, sampleSize));
	EXGskew = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2);
	Beta = w' * EXGskew;
	w = w - myy * (EXGskew - Beta*w)/(-Beta);
	
       otherwise
	error('Code for desired nonlinearity not found!');
      end
      
      % Normalize the new w.
      w = w / norm(w);
      i = i + 1;
    end
    round = round + 1;
  end
  if b_verbose, fprintf('Done.\n'); end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Also plot the ones that may not have been plotted.
  if (usedDisplay > 0) & (rem(round-1, displayInterval) ~= 0)
    switch usedDisplay
     case 1
      % There was and may still be other displaymodes...
      % 1D signals
      temp = X'*B;
      icaplot(plottype,temp(:,1:numOfIC)');
      drawnow;
     case 2
      % ... and now there are :-)
      % 1D basis
      icaplot(plottype,A');
      %drawnow;
     case 3
      % ... and now there are :-)
      % 1D filters
      icaplot(plottype,W);
      %drawnow;
     otherwise
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the end let's check the data for some security
if imag(A) <> 0
  if b_verbose, fprintf('Warning: removing the imaginary part from the result.\n'); end
  A = real(A);
  W = real(W);
end

endfunction

                                                                                                                                                                                                                                                                                getSamples.m                                                                                        0100644 0017072 0006200 00000000146 07450203724 011777  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function Samples = getSamples(max, percentage)
Samples = find(rand(1, max) < percentage);
endfunction
                                                                                                                                                                                                                                                                                                                                                                                                                          icaplot.m                                                                                           0100644 0017072 0006200 00000026243 07450203724 011334  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function icaplot(mode, ...);

%ICAPLOT - plot signals in various ways
%
% ICAPLOT is mainly for plottinf and comparing the mixed signals and
% separated ica-signals.
%
% ICAPLOT has many different modes. The first parameter of the function
% defines the mode. Other parameters and their order depends on the
% mode. The explanation for the more common parameters is in the end.
%
% Classic
%     icaplot('classic', s1, n1, range, xrange, titlestr)
%
%     Plots the signals in the same manner as the FASTICA and FASTICAG
%     programs do. All the signals are plotted in their own axis.
%
% Complot
%     icaplot('complot', s1, n1, range, xrange, titlestr)
%
%     The signals are plotted on the same axis. This is good for
%     visualization of the shape of the signals. The scale of the signals 
%     has been altered so that they all fit nicely.
%
% Histogram
%     icaplot('histogram', s1, n1, range, bins, style)
%     
%     The histogram of the signals is plotted. The number of bins can be
%     specified with 'bins'-parameter. The style for the histograms can
%     be either 'bar' (default) of 'line'.
%
% Scatter
%     icaplot('scatter', s1, n1, s2, n2, range, titlestr, s1label,
%     s2label, markerstr)
%
%     A scatterplot is plotted so that the signal 1 is the 'X'-variable
%     and the signal 2 is the 'Y'-variable. The 'markerstr' can be used
%     to specify the maker used in the plot. The format for 'markerstr'
%     is the same as for Matlab's PLOT. 
%
% Compare
%     icaplot('compare', s1, n1, s2, n2, range, xrange, titlestr,
%     s1label, s2label)
%
%     This is for comparing two signals. The main used in this context
%     would probably be to see how well the separated ICA-signals explain 
%     the observed mixed signals. The s2 signals are first scaled with
%     REGRESS function.
%
% Compare - Sum
%     icaplot('sum', s1, n1, s2, n2, range, xrange, titlestr, s1label,
%     s2label)
%
%     The same as Compare, but this time the signals in s2 (specified by
%     n2) are summed together.
%
% Compare - Sumerror
%     icaplot('sumerror', s1, n1, s2, n2, range, xrange, titlestr,
%     s1label, s2label)
%     
%     The same as Compare - Sum, but also the 'error' between the signal
%     1 and the summed IC's is plotted.
%
%
% More common parameters
%     The signals to be plotted are in matrices s1 and s2. The n1 and n2
%     are used to tell the index of the signal or signals to be plotted
%     from s1 or s2. If n1 or n2 has a value of 0, then all the signals
%     from corresponding matrix will be plotted. The values for n1 and n2 
%     can also be vectors (like: [1 3 4]) In some casee if there are more
%     than 1 signal to be plotted from s1 or s2 then the plot will
%     contain as many subplots as are needed. 
%
%     The range of the signals to be plotted can be limited with
%     'range'-parameter. It's value is a vector ( 10000:15000 ). If range 
%     is 0, then the whole range will be plotted.
%
%     The 'xrange' is used to specify only the labels used on the
%     x-axis. The value of 'xrange' is a vector containing the x-values
%     for the plots or [start end] for begin and end of the range
%     ( 10000:15000 or [10 15] ). If xrange is 0, then value of range
%     will be used for x-labels.
%
%     You can give a title for the plot with 'titlestr'. Also the
%     's1label' and 's2label' are used to give more meaningfull label for 
%     the signals.
%
%     Lastly, you can omit some of the arguments from the and. You will
%     have to give values for the signal matrices (s1, s2) and the
%     indexes (n1, n2)

% 7.8.1998

% Added by DR:
% Octave complains if we run icaplot.m as a script file.  So
% call the rest of the functions needed in the following script:

icaplotfunctions;

str_param = lower(mode);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'dispsig' is to replace the old DISPSIG
% '' & 'classic' are just another names - '' quite short one :-)

if strcmp(str_param, '') | strcmp(str_param, 'classic') | strcmp(str_param,'dispsig')
  
  % icaplot(mode, s1, n1, range, xrange, titlestr)
  va_start();
  if nargin-1 < 1, error('Not enough arguments.'); endif
  s1 = va_arg();
  if nargin-1 < 2, n1 = 0; else n1 = va_arg(); endif
  if nargin-1 < 3, range = 0;else range = va_arg(); endif
  if nargin-1 < 4, xrange = 0;else xrange = va_arg(); endif
  if nargin-1 < 5, titlestr = ''; else titlestr = va_arg(); endif
  range=chkrange(range, s1);
  xrange=chkxrange(xrange, range);
  n1=chkn(n1, s1);
  clg;
  
  numSignals = size(n1, 2);
  for i = 1:numSignals,
    subplot(numSignals, 1, i);
    % Added by DR
    clg;
    % "if" statement added by DR to prevent from trying to plot empty matrices
  if (!all(s1(n1(i),range) == zeros(size(range)))),
    plot(xrange, s1(n1(i), range)); 
  endif
  end
  subplot(numSignals,1, 1);
  if (~isempty(titlestr))
    title(titlestr);
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(str_param, 'complot'),
  % icaplot(mode, s1, n1, range, xrange, titlestr)
  va_start();
  if nargin-1 < 1, error('Not enough arguments.'); end
  s1 = remmean(va_arg());
  if nargin-1 < 2, n1 = 0;else n1 = va_arg(); end
  if nargin-1 < 3, range = 0;else range = va_arg(); end
  if nargin-1 < 4, xrange = 0;else xrange = va_arg(); end
  if nargin-1 < 5, titlestr = '';else titlestr = va_arg(); end
  
  range=chkrange(range, s1);
  xrange=chkxrange(xrange, range);
  n1=chkn(n1, s1);
  
  for i = 1:size(n1, 2)
    S1(i, :) = s1(n1(i), range);
  end
  
  alpha = mean(max(S1')-min(S1'));
  for i = 1:size(n1,2)
    S2(i,:) = S1(i,:) - alpha*(i-1)*ones(size(S1(1,:)));
  end
  
  plot(xrange, S2');
  axis([min(xrange) max(xrange) min(min(S2)) max(max(S2)) ]);
  
  set(gca,'YTick',(-size(S1,1)+1)*alpha:alpha:0);
  set(gca,'YTicklabel',fliplr(n1));
  
  if (~isempty(titlestr))
    title(titlestr);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(str_param, 'histogram'),
  % icaplot(mode, s1, n1, range, bins, style)
  if nargin-1 < 1, error('Not enough arguments.'); end
  s1 = va_arg();
  if nargin-1 < 2, n1 = 0;else n1 = va_arg(); end
  if nargin-1 < 3, range = 0;else range = va_arg(); end
  if nargin-1 < 4, bins = 100;else bins = va_arg(); end
  if nargin-1 < 5, style = 'bar';else style = va_arg(); end

  range = chkrange(range, s1);
  n1 = chkn(n1, s1);
  % Added by DR
  % Make the number of bins a factor of ten different from 
  % the number of data points:

  numSignals = size(n1, 2);
  rows = floor(sqrt(numSignals));
  columns = ceil(sqrt(numSignals));
  while (rows * columns < numSignals)
    columns = columns + 1;
  end
  
  str_param = lower(style);
  if strcmp(str_param, 'bar'),
    for i = 1:numSignals,
      subplot(rows, columns, i);
      % Added by DR:
      clg;
    % "if" statement added by DR to prevent from trying to plot empty matrices
    if (!all(s1(n1(i),range) == zeros(size(range)))),
       hist(s1(n1(i), range), bins);
    endif  
      title(int2str(n1(i)));
      %drawnow;
    endfor

  elseif strcmp(str_param,''),
    for i = 1:numSignals,
      subplot(rows, columns, i);
      % Added by DR:
      clg;
 % "if" statement added by DR to prevent from trying to plot empty matrices
    if (!all(s1(n1(i),range) == zeros(size(range)))),
      hist(s1(n1(i), range), bins);
    endif
      title(int2str(n1(i)));
      %drawnow;
    end

  elseif strcmp(str_param, 'line'),
    for i = 1:numSignals,
      subplot(rows, columns, i);
    if (!all(s1(n1(i),range) == zeros(size(range)))),
      [Y, X]=hist(s1(n1(i), range), bins);
    endif
      % Added by DR:
      clg;
      %plot(X, Y);
      title(int2str(n1(i)));
      %drawnow;
    end
  else
    fprintf('Unknown style.\n')
  endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(str_param, 'scatter'),
  % icaplot(mode, s1, n1, s2, n2, range, titlestr, xlabelstr, ylabelstr, markerstr)
  if nargin-1 < 4, error('Not enough arguments.'); end
  s1 = va_arg();
  n1 = va_arg();
  s2 = va_arg();
  n2 = va_arg();
  if nargin-1 < 5, range = 0;else range = va_arg(); end
  if nargin-1 < 6, titlestr = '';else titlestr = va_arg(); end
  if nargin-1 < 7, xlabelstr = 'Signal 1';else xlabelstr = va_arg(); end
  if nargin-1 < 8, ylabelstr = 'Signal 2';else ylabelstr = va_arg(); end
  if nargin-1 < 9, markerstr = '.';else markerstr = va_arg(); end

  range = chkrange(range, s1);
  n1 = chkn(n1, s1);
  n2 = chkn(n2, s2);
  
  rows = size(n1, 2);
  columns = size(n2, 2);
  for r = 1:rows
    for c = 1:columns
      subplot(rows, columns, (r-1)*columns + c);
      % Added by DR:
      clg;
      plot(s1(n1(r), range),s2(n2(c), range),markerstr);
      if (~isempty(titlestr))
	title(titlestr);
      end
      if (rows*columns == 1)
	xlabel(xlabelstr);
	ylabel(ylabelstr);
      else 
	xlabel([xlabelstr ' (' int2str(n1(r)) ')']);
	ylabel([ylabelstr ' (' int2str(n2(c)) ')']);
      end
      %drawnow;
    end
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(str_param, 'compare') | strcmp(str_param, 'sum') | strcmp(str_param,'sumerror'),
  % icaplot(mode, s1, n1, s2, n2, range, xrange, titlestr, s1label, s2label)
  va_start();
  if nargin-1 < 4, error('Not enough arguments.'); end
  s1 = va_arg();
  n1 = va_arg();
  s2 = va_arg();
  n2 = va_arg();
  if nargin-1 < 5, range = 0;else range = va_arg(); end
  if nargin-1 < 6, xrange = 0;else xrange = va_arg(); end
  if nargin-1 < 7, titlestr = '';else titlestr = va_arg(); end
  if nargin-1 < 8, s1label = 'Mix';else s1label = va_arg(); end
  if nargin-1 < 9, s2label = 'IC';else s2label = va_arg(); end

  range = chkrange(range, s1);
  xrange = chkxrange(xrange, range);
  n1 = chkn(n1, s1);
  n2 = chkn(n2, s2);

  numSignals = size(n1, 2);
  if (numSignals > 1)
    externalLegend = 1;
  else
    externalLegend = 0;
  end
  
  rows = floor(sqrt(numSignals+externalLegend));
  columns = ceil(sqrt(numSignals+externalLegend));
  while (rows * columns < (numSignals+externalLegend))
    columns = columns + 1;
  end
  
  clf;
  
  for j = 1:numSignals
    subplot(rows, columns, j);
    str_param = lower(mode);
    if strcmp(str_param, 'compare')
      plotcompare(s1, n1(j), s2,n2, range, xrange);
      [legendtext,legendstyle]=legendcompare(n1(j),n2,s1label,s2label,externalLegend);
    elseif strcmp(str_param, 'sum')
      plotsum(s1, n1(j), s2,n2, range, xrange);
      [legendtext,legendstyle]=legendsum(n1(j),n2,s1label,s2label,externalLegend);
    elseif strcmp(str_param, 'sumerror')
      plotsumerror(s1, n1(j), s2,n2, range, xrange);
      [legendtext,legendstyle]=legendsumerror(n1(j),n2,s1label,s2label,externalLegend);
    endif
    
    if externalLegend
      title([titlestr ' (' s1label  ' ' int2str(n1(j)) ')']);
    else
      legend(char(legendtext));
      if (~isempty(titlestr))
	title(titlestr);
      end
    end
  end
  
  if (externalLegend)
    subplot(rows, columns, numSignals+1);
    legendsize = size(legendtext, 2);
    hold on;
    for i=1:legendsize
      plot([0 1],[legendsize-i legendsize-i], char(legendstyle(i)));
      text(1.5, legendsize-i, char(legendtext(i)));
    end
    hold off;
    axis([0 6 -1 legendsize]);
    %axis off;
  end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
endif

endfunction


                                                                                                                                                                                                                                                                                                                                                             icaplotfunctions.m                                                                                  0100644 0017072 0006200 00000006012 07450203724 013255  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    % This is a script file called from "icaplot.m"
1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotcompare(s1, n1, s2, n2, range, xrange);
  style=getStyles;
  K = regress(s1(n1,:)',s2');
  plot(xrange, s1(n1,range), char(style(1)));
  hold on
  for i=1:size(n2,2)
    plotstyle=char(style(i+1));
    plot(xrange, K(n2(i))*s2(n2(i),range), plotstyle);
  end
  hold off
endfunction

function [legendText, legendStyle]=legendcompare(n1, n2, s1l, s2l, externalLegend);
  style=getStyles;
  if (externalLegend)
    legendText(1)=[s1l ' (see the titles)'];
  else
    legendText(1)=[s1l ' ', int2str(n1)];
  end
  legendStyle(1)=style(1);
  for i=1:size(n2, 2)
    legendText(i+1) = [s2l ' ' int2str(n2(i))];
    legendStyle(i+1) = style(i+1);
  end
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotsum(s1, n1, s2, n2, range, xrange);
  K = diag(regress(s1(n1,:)',s2'));
  sigsum = sum(K(:,n2)*s2(n2,:));
  plot(xrange, s1(n1, range),'k-', ...
       xrange, sigsum(range), 'b-');
endfunction

function [legendText, legendStyle]=legendsum(n1, n2, s1l, s2l, externalLegend);
  if (externalLegend)
    legendText(1)=[s1l ' (see the titles)'];
  else
    legendText(1)=[s1l ' ', int2str(n1)];
  end
  legendText(2)=['Sum of ' s2l ': ', int2str(n2)];
  legendStyle=['k-';'b-'];
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotsumerror(s1, n1, s2, n2, range, xrange);
  K = diag(regress(s1(n1,:)',s2'));
  sigsum = sum(K(:,n2)*s2(n2,:));
  plot(xrange, s1(n1, range),'k-', ...
       xrange, sigsum(range), 'b-', ...
       xrange, s1(n1, range)-sigsum(range), 'r-');
endfunction

function [legendText, legendStyle]=legendsumerror(n1, n2, s1l, s2l, externalLegend);
  if (externalLegend)
    legendText(1)=[s1l ' (see the titles)'];
  else
    legendText(1)=[s1l ' ', int2str(n1)];
  end
  legendText(2)=['Sum of ' s2l ': ', int2str(n2)];
  legendText(3)='"Error"';
  legendStyle=['k-';'b-';'r-'];
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function style=getStyles;
  color = ['k','r','g','b','m','c','y'];
  line = ['-',':','-.','--'];
  for i = 0:size(line,2)-1
    for j = 1:size(color, 2)
      style(j + i*size(color, 2)) = strcat(color(j), line(i+1));
    end
  end
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function range=chkrange(r, s)
  if r == 0
    range = 1:size(s, 2);
  else
    range = r;
  end
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xrange=chkxrange(xr,r);
  if xr == 0
    xrange = r;
  elseif size(xr, 2) == 2
    xrange = xr(1):(xr(2)-xr(1))/(size(r,2)-1):xr(2);
  elseif size(xr, 2)~=size(r, 2)
    error('Xrange and range have different sizes.');
  else
    xrange = xr;
  end
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n=chkn(n,s)
  if n == 0
    n = 1:size(s, 1);
  end
endfunction
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      pcamat.m                                                                                            0100644 0017072 0006200 00000026524 07450203724 011150  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function [E, D] = pcamat(vectors, firstEig, lastEig, s_interactive, ...
    s_verbose);
%PCAMAT - Calculates the pca for data
%
% [E, D] = pcamat(vectors, firstEig, lastEig, ... 
%                 interactive, verbose);
%
% Calculates the PCA matrices for given data (row) vectors. Returns
% the eigenvector (E) and diagonal eigenvalue (D) matrices containing the
% selected subspaces. Dimensionality reduction is controlled with
% the parameters 'firstEig' and 'lastEig' - but it can also be done
% interactively by setting parameter 'interactive' to 'on' or 'gui'.
%
% ARGUMENTS
%
% vectors       Data in row vectors.
% firstEig      Index of the largest eigenvalue to keep.
%               Default is 1.
% lastEig       Index of the smallest eigenvalue to keep.
%               Default is equal to dimension of vectors.
% interactive   Specify eigenvalues to keep interactively. Note that if
%               you set 'interactive' to 'on' or 'gui' then the values
%               for 'firstEig' and 'lastEig' will be ignored, but they
%               still have to be entered. If the value is 'gui' then the
%               same graphical user interface as in FASTICAG will be
%               used. Default is 'off'.
% verbose       Default is 'on'.
%
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%
% Note 
%       The eigenvalues and eigenvectors returned by PCAMAT are not sorted.
%
% This function is needed by FASTICA and FASTICAG

% For historical reasons this version does not sort the eigenvalues or
% the eigen vectors in any ways. Therefore neither does the FASTICA or
% FASTICAG. Generally it seams that the components returned from
% whitening is almost in reversed order. (That means, they usually are,
% but sometime they are not - depends on the EIG-command of matlab.)

% 16.6.2000
% Hugo G?vert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values:
if nargin < 5, s_verbose = 'on'; end
if nargin < 4, s_interactive = 'off'; end
if nargin < 3, lastEig = size(vectors, 1); end
if nargin < 2, firstEig = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the optional parameters;
if strcmp(lower(s_verbose), 'on')
    b_verbose = 1;
  elseif strcmp(lower(s_verbose), 'off')
    b_verbose = 0;
  else
    error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

if strcmp(lower(s_interactive), 'on')
    b_interactive = 1;
  elseif strcmp(lower(s_interactive), 'off')
    b_interactive = 0;
  elseif strcmp(lower(s_interactive), 'gui')
    b_interactive = 2;
  else
    error(sprintf('Illegal value [ %s ] for parameter: ''interactive''\n', ...
          s_interactive));
end

oldDimension = size (vectors, 1);
if ~(b_interactive)
  if lastEig < 1 | lastEig > oldDimension
    error(sprintf('Illegal value [ %d ] for parameter: ''lastEig''\n', lastEig));
  end
  if firstEig < 1 | firstEig > lastEig
    error(sprintf('Illegal value [ %d ] for parameter: ''firstEig''\n', firstEig));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate PCA

% Calculate the covariance matrix.
if b_verbose, fprintf ('Calculating covariance...\n'); end
covarianceMatrix = cov(vectors')*(size(vectors',1)-1)/size(vectors',1);

maxLastEig = rank(covarianceMatrix, 1e-9);

% Calculate the eigenvalues and eigenvectors of covariance matrix.
[E, D] = eig(covarianceMatrix);

% Sort the eigenvalues - decending.
eigenvalues = flipud(sort(diag(D)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive part - command-line
if b_interactive == 1

  % Show the eigenvalues to the user
  hndl_win=figure;
  bar(eigenvalues);
  title('Eigenvalues');

  % ask the range from the user...
  % ... and keep on asking until the range is valid :-)
  areValuesOK=0;
  while areValuesOK == 0
    firstEig = input('The index of the largest eigenvalue to keep? (1) ');
    lastEig = input(['The index of the smallest eigenvalue to keep? (' ...
                    int2str(oldDimension) ') ']);
    % Check the new values...
    % if they are empty then use default values
    if isempty(firstEig), firstEig = 1;end
    if isempty(lastEig), lastEig = oldDimension;end
    % Check that the entered values are within the range
    areValuesOK = 1;
    if lastEig < 1 | lastEig > oldDimension
      fprintf('Illegal number for the last eigenvalue.\n');
      areValuesOK = 0;
    end
    if firstEig < 1 | firstEig > lastEig
      fprintf('Illegal number for the first eigenvalue.\n');
      areValuesOK = 0;
    end
  end
  % close the window
  close(hndl_win);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive part - GUI
if b_interactive == 2

  % Show the eigenvalues to the user
  hndl_win = figure('Color',[0.8 0.8 0.8], ...
    'PaperType','a4letter', ...
    'Units', 'normalized', ...
    'Name', 'FastICA: Reduce dimension', ...
    'NumberTitle','off', ...
    'Tag', 'f_eig');
  h_frame = uicontrol('Parent', hndl_win, ...
    'BackgroundColor',[0.701961 0.701961 0.701961], ...
    'Units', 'normalized', ...
    'Position',[0.13 0.05 0.775 0.17], ...
    'Style','frame', ...
    'Tag','f_frame');

b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.142415 0.0949436 0.712077 0.108507], ...
	'String','Give the indices of the largest and smallest eigenvalues of the covariance matrix to be included in the reduced data.', ...
	'Style','text', ...
	'Tag','StaticText1');
e_first = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'f=round(str2num(get(gcbo, ''String'')));' ...
          'if (f < 1), f=1; end;' ...
          'l=str2num(get(findobj(''Tag'',''e_last''), ''String''));' ...
          'if (f > l), f=l; end;' ...
          'set(gcbo, ''String'', int2str(f));' ...
          ], ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'Position',[0.284831 0.0678168 0.12207 0.0542535], ...
	'Style','edit', ...
        'String', '1', ...
	'Tag','e_first');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.142415 0.0678168 0.12207 0.0542535], ...
	'String','Range from', ...
	'Style','text', ...
	'Tag','StaticText2');
e_last = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'l=round(str2num(get(gcbo, ''String'')));' ...
          'lmax = get(gcbo, ''UserData'');' ...
          'if (l > lmax), l=lmax; fprintf([''The selected value was too large, or the selected eigenvalues were close to zero\n'']); end;' ...
          'f=str2num(get(findobj(''Tag'',''e_first''), ''String''));' ...
          'if (l < f), l=f; end;' ...
          'set(gcbo, ''String'', int2str(l));' ...
          ], ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'Position',[0.467936 0.0678168 0.12207 0.0542535], ...
	'Style','edit', ...
        'String', int2str(maxLastEig), ...
        'UserData', maxLastEig, ...
	'Tag','e_last');
% in the first version oldDimension was used instead of 
% maxLastEig, but since the program would automatically
% drop the eigenvalues afte maxLastEig...
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.427246 0.0678168 0.0406901 0.0542535], ...
	'String','to', ...
	'Style','text', ...
	'Tag','StaticText3');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback','uiresume(gcbf)', ...
	'Position',[0.630697 0.0678168 0.12207 0.0542535], ...
	'String','OK', ...
	'Tag','Pushbutton1');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'gui_help(''pcamat'');' ...
          ], ...
	'Position',[0.767008 0.0678168 0.12207 0.0542535], ...
	'String','Help', ...
	'Tag','Pushbutton2');

  h_axes = axes('Position' ,[0.13 0.3 0.775 0.6]);
  set(hndl_win, 'currentaxes',h_axes);
  bar(eigenvalues);
  title('Eigenvalues');

  uiwait(hndl_win);
  firstEig = str2num(get(e_first, 'String'));
  lastEig = str2num(get(e_last, 'String'));

  % close the window
  close(hndl_win);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See if the user has reduced the dimension enought

if lastEig > maxLastEig
  lastEig = maxLastEig;
  if b_verbose
    fprintf('Dimension reduced to %d due to the singularity of covariance matrix\n',...
           lastEig-firstEig+1);
  end
else
  % Reduce the dimensionality of the problem.
  if b_verbose
    if oldDimension == (lastEig - firstEig + 1)
      fprintf ('Dimension not reduced.\n');
    else
      fprintf ('Reducing dimension...\n');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drop the smaller eigenvalues
if lastEig < oldDimension
  lowerLimitValue = (eigenvalues(lastEig) + eigenvalues(lastEig + 1)) / 2;
else
  lowerLimitValue = eigenvalues(oldDimension) - 1;
end

lowerColumns = diag(D) > lowerLimitValue;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drop the larger eigenvalues
if firstEig > 1
  higherLimitValue = (eigenvalues(firstEig - 1) + eigenvalues(firstEig)) / 2;
else
  higherLimitValue = eigenvalues(1) + 1;
end
higherColumns = diag(D) < higherLimitValue;

% Combine the results from above
selectedColumns = lowerColumns & higherColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print some info for the user
if b_verbose
  fprintf ('Selected [ %d ] dimensions.\n', sum (selectedColumns));
end
if sum (selectedColumns) ~= (lastEig - firstEig + 1),
  error ('Selected a wrong number of dimensions.');
end

if b_verbose
  fprintf ('Smallest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(lastEig));
  fprintf ('Largest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(firstEig));
  fprintf ('Sum of removed eigenvalues [ %g ]\n', sum(diag(D) .* ...
    (~selectedColumns)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the colums which correspond to the desired range
% of eigenvalues.
E = selcol (E, selectedColumns);
D = selcol (selcol (D, selectedColumns)', selectedColumns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some more information
if b_verbose
  sumAll=sum(eigenvalues);
  sumUsed=sum(diag(D));
  retained = (sumUsed / sumAll) * 100;
  fprintf('[ %g ] %% of (non-zero) eigenvalues retained.\n', retained);
end
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Moved to "selcol.m"

%function newMatrix = selcol(oldMatrix, maskVector);

% newMatrix = selcol(oldMatrix, maskVector);

% Selects the columns of the matrix that marked by one in the given vector.
% The maskVector is a column vector.

% 15.3.1998

%if size(maskVector, 1) ~= size(oldMatrix, 2),
%  error ('The mask vector and matrix are of uncompatible size.');
%end

%numTaken = 0;

%for i = 1 : size (maskVector, 1),
%  if maskVector(i, 1) == 1,
%    takingMask(1, numTaken + 1) = i;
%    numTaken = numTaken + 1;
%  end
%end

%newMatrix = oldMatrix(:, takingMask);
%endfunction
                                                                                                                                                                            remmean.m                                                                                           0100644 0017072 0006200 00000000654 07450203724 011323  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function [newVectors, meanValue] = remmean(vectors);
%REMMEAN - remove the mean from vectors
%
% [newVectors, meanValue] = remmean(vectors);
%
% Removes the mean of row vectors.
% Returns the new vectors and the mean.
%
% This function is needed by FASTICA and FASTICAG

% 24.8.1998
% Hugo G?vert

newVectors = zeros (size (vectors));
meanValue = mean (vectors')';
newVectors = vectors - meanValue * ones (1,size (vectors, 2));
                                                                                    selcol.m                                                                                            0100644 0017072 0006200 00000001043 07450203724 011151  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function newMatrix = selcol(oldMatrix, maskVector);

% newMatrix = selcol(oldMatrix, maskVector);

% Selects the columns of the matrix that marked by one in the given vector.
% The maskVector is a column vector.

% 15.3.1998

if size(maskVector, 1) ~= size(oldMatrix, 2),
  error ('The mask vector and matrix are of uncompatible size.');
end

numTaken = 0;

for i = 1 : size (maskVector, 1),
  if maskVector(i, 1) == 1,
    takingMask(1, numTaken + 1) = i;
    numTaken = numTaken + 1;
  end
end

newMatrix = oldMatrix(:, takingMask);
endfunction
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             whitenv.m                                                                                           0100644 0017072 0006200 00000004064 07450203724 011362  0                                                                                                    ustar   dryan                           cdf                                                                                                                                                                                                                    function [newVectors, whiteningMatrix, dewhiteningMatrix] = whitenv ...
    (vectors, E, D, s_verbose);
%WHITENV - Whitenv vectors.
%
% [newVectors, whiteningMatrix, dewhiteningMatrix] = ...
%                               whitenv(vectors, E, D, verbose);
%
% Whitens the data (row vectors) and reduces dimension. Returns
% the whitened vectors (row vectors), whitening and dewhitening matrices.
%
% ARGUMENTS
%
% vectors       Data in row vectors.
% E             Eigenvector matrix from function 'pcamat'
% D             Diagonal eigenvalue matrix from function 'pcamat'
% verbose       Optional. Default is 'on'
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%       [nv, wm, dwm] = whitenv(vectors, E, D);
%
%
% This function is needed by FASTICA and FASTICAG
%
%   See also PCAMAT

% 24.8.1998
% Hugo G?vert

% ========================================================
% Default value for 'verbose'
if nargin < 4, s_verbose = 'on'; end

% Check the optional parameter verbose;
switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

% ========================================================
% Calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously).
whiteningMatrix = inv (sqrt (D)) * E';
dewhiteningMatrix = E * sqrt (D);

% Project to the eigenvectors of the covariance matrix.
% Whiten the samples and reduce dimension simultaneously.
if b_verbose, fprintf ('Whitening...\n'); end
newVectors =  whiteningMatrix * vectors;

% ========================================================
% Just some security...
if max (max (imag (newVectors))) ~= 0,
  error ('Whitened vectors have imaginary values.');
end

% Print some information to user
if b_verbose
  fprintf ('Check: covariance differs from identity by [ %g ].\n', ...
% Added by DR for "Octave" compatibility:
% Normalize the covariance matrix by N, not by N-1:
max (max (abs (cov(newVectors')*((size(newVectors',1)-1)/size(newVectors',1)) - eye (size (newVectors, 1))))));
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            