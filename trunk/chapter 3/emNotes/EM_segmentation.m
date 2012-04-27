function [mask,mu,v,p]=EM_segmentation(image,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Expectation Maximization image segmentation
%
%   Input:
%          image: grey color image
%          M: Number of classes
%   Output:
%          mask: clasification image mask
%          mu: vector of class means 
%          v: vector of class variances
%          p: vector of class proportions (weights)   
%
%   Example: [mask,mu,v,p]=EM_segmentation(image,3);
%
%   Author: Prof. Jose Vicente Manjon Herrera
%    Email: jmanjon@fis.upv.es
%     Date: 02-05-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check image
image        = double(image);
copy         = image;          
min_orig     = min(image(:));     
image        = image - min_orig + 1;  % Maps image so lowest gray level goes to one.
max_scal_img = max(image(:));            


% Create image histogram

Px = hist(image(:) , ceil(max_scal_img)+1 ); % Calc. histogram with max. num. of bins as max. gray level.
Px = conv(Px , [1,2,3,2,1]);          % Smoothing the histogram with linear Bspline.
%Px = conv(Px , [1,3,6,7,6,3,1]);
Px = Px(3:(length(Px)-2));            % Clipping out the filter support regions. 
%Px = Px(4:(length(Px)-3));            % Clipping out the filter support regions. 
Px = Px/sum(Px);                      % Transforming histogram to discrete pdf.

x  = find(Px)';  % Taking only the gray levels that are present in the image (histogram bins).
Px = Px(x)';     % Taking the normalized counts for each of the gray level bins.

% Initial parameters

mu = (1:M)*max_scal_img/(M+1); % Mixture Component Initial Means.
v  = ones(1,M)*max_scal_img;   % Mixture Component Initial Variances.
p  = ones(1,M)*1/M;            % Mixture Component Initial Weights (must add to one).


% Expectation-Maximization process

%figure;
%plot(x,Px,'b','LineWidth',1.5);

%hold on
	%plot(x,gauss_pdf(mu,v,p,x),'g:','LineWidth',1.5) ;
%plot(x,sum(gauss_pdf(mu,v,p,x),2),'-.r','LineWidth',1.5);
%drawnow
% 
% figure;

sml = mean(diff(x))/1000;
k = 0;
while(1)
    k = k + 1;
    % Expectation Part
    prb    = gauss_pdf(mu,v,p,x);  % alpha_l * p_l(x_i | Theta_l_g) where l=[1..M] and i=[1..#gray level bins]
    scal   = sum(prb,2) + eps;     % Sum_from_1_to_#samples[ alpha_l * p_l(x_i | Theta_l_g) ] o
    loglik = sum(Px.*log(scal));   % Incomplete-data likelihood
    
    % Maximization Part
    for j=1:M                      % Repeat for each mixture component (l).
        pp   = Px.*prb(:,j)./scal;  
        p(j) = sum(pp);
        mu(j)= sum(x.*pp)/p(j);
        vr   = (x-mu(j));
        v(j) = sum(vr.*vr.*pp)/p(j)+sml;
    end
    p = p + 1e-3;
    p = p/sum(p);

    % Exit Condition (log-likelihood comparison)
    prb     = gauss_pdf(mu,v,p,x);
    scal    = sum(prb,2)+eps;
    nloglik = sum(Px.*log(scal));
    
     %clf
     %figure
	plot (x,Px);
	print ("plot.png", '-dpng');

     %plot(x,Px,'b','LineWidth',1.5);
     %hold on
     %plot(x,prb,'g:','LineWidth',1.5)
	plot (x, prb);
	print ("probability.png", '-dpng');
     %plot(x,sum(prb,2),'-.r','LineWidth',1.5)
		plot (x, sum(prb,2));
		print ("probability-2.png", '-dpng');
     %drawnow
    
    if((nloglik-loglik)<0.0001) break; end;        
end

% Calculates mask

mu       = mu + min_orig - 1;        % Maps the estimated means into the original image space
imsize   = size(copy);
mask_aux = zeros([imsize,M]);

for j=1:M
    mask_aux(:,:,j) = (ones(imsize)*p(j)) .* normpdf( copy , ones(imsize)*mu(j) , sqrt(ones(imsize)*v(j)) );
end
[max_val , mask] = max(mask_aux , [] ,3);



%%%%%%%%%%%%%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%

function y = gauss_pdf(m,v,g,x)
x=x(:);
m=m(:);
v=v(:);
g=g(:);
for i=1:size(m,1)
   d = x-m(i);
   amp = g(i)/sqrt(2*pi*v(i));
   y(:,i) = amp*exp(-0.5 * (d.*d)/v(i));
end


