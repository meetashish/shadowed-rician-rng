function [f,F] = EstimateDistribution(X,x)
% This function implements estimation of CDF and PDF of one dimensional 
% random variables.
%
% INPUTS:
%           X = vector specifying random variables
%           x = vector specifying points for which CDF and PDF has to be
%               evaluated
% OUTPUTS:
%           f = vector specifying estimated PDF of random variable X for
%               points.
%           F = vector specifying estimated CDF of random variable X for
%               points.
%
% USAGE EXAMPLES:
% %% Generate N Standard Normally Distributed Random Variable
% N = 1000000;
% X = randn(N,1);
% %% Points for which CDF and PDF are to be evaluated
% x = linspace(-10,10,1000);
% %% Estimate PDF and CDF
% [f,F] = EstimateDistribution(X,x);
% %% Plot Results
% figure(1);
% plot(x,f,x,F);
% xlabel('x');
% ylabel('Simulated PDF & CDF');
% str1 = strcat('PDF;','Area = ',num2str(trapz(x,f)));
% legend(str1,'CDF','Location','northwest');
% 
% 
% %% Generate N Gaussianly Distributed Random Variable with specific mean and
% %% Standard Deviation
% N = 1000000;
% mu = -1;
% sigma = 5;
% X = mu + sigma*randn(N,1);
% %% Points for which CDF and PDF are to be evaluated
% x = linspace(-10,10,1000);
% %% Theoretical PDF and CDF
% fx = (1/sqrt((2*pi*sigma*sigma)))*exp(-(((x - mu).^2)/(2*sigma*sigma)));
% Fx = 0.5*(1 + erf((x - mu)/(sqrt(2*sigma*sigma))));
% %% Estimate PDF and CDF
% [f,F] = EstimateDistribution(X,x);
% %% Plot Results
% figure(2);
% plot(x,f,x,fx,x,F,x,Fx);
% xlabel('x');
% ylabel('PDF & CDF');
% str1 = strcat('Simulated PDF;','Area = ',num2str(trapz(x,f)));
% str2 = strcat('Theoretical PDF;','Area = ',num2str(trapz(x,fx)));
% legend(str1,str2,'Simulated CDF','Theoretical CDF','Location','northwest');
% 
% 
% %% Generate N Uniformaly Distributed Random Variable
% N = 1000000;
% X = rand(N,1);
% %% Points for which CDF and PDF are to be evaluated
% x = linspace(-10,10,1000);
% %% Estimate PDF and CDF
% [f,F] = EstimateDistribution(X,x);
% %% Plot Results
% figure(3);
% plot(x,f,x,F);
% xlabel('x');
% ylabel('Simulated PDF & CDF');
% str1 = strcat('PDF;','Area = ',num2str(trapz(x,f)));
% legend(str1,'CDF','Location','northwest');
%
% REFERENCES:
% Athanasios Papoulis, S. Unnikrishna Pillai, Probability, Random Variables
% and Stochastic Processes, 4e
% Peyton Z. Peebles Jr., Probability, Random Variables, And Random Signal 
% Principles, 2e
% Saeed Ghahramani, Fundamentals of Probability, with Stochastic Processes,
% 3e
%
% SEE ALSO:
% interp1, smooth
%
% AUTHOR:
% Ashish (Meet) Meshram
% meetashish85@gmail.com; mt1402102002@iiti.ac.in

% Checking Input Arguments
if nargin<2||isempty(x), x = linspace(-10,10,1000);end
if nargin<2||isempty(X)
    error('Missing Input Arguments: Please specify vector random variables');
end

% Impelementation Starts Here
f = zeros(1,length(x)); % Preallocation of memory space
F = f;                  % Preallocation of memory space
h = 0.000000001;        % Small value closer to zero for evaluating
                        % numerical differentiation.

% Estimating CDF by its definition
for m = 1:length(x)
    p = 0;              % True Probability
    q = 0;              % False Probability
    for n = 1:length(X)
        if X(n)<=x(m)   % Definition of CDF
            p = p + 1;
        else
            q = q + 1;
        end
    end
    F(m) = p/(p + q);   % Calulating Probability
end

% Estimating PDF by differentiation of CDF
for k = 1:length(x)
    fxph = interp1(x,F,x(k) + h,'spline');  % Interpolating value of F(x+h)
    fxmh = interp1(x,F,x(k) - h,'spline');  % Interpolating value of F(x-h)
    f(k) = (fxph - fxmh)/(2*h);             % Two-point formula to compute
end                                         % Numerical differentiation
f = smooth(f);                              % Smoothing at last