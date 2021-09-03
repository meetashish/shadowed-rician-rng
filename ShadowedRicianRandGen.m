function X = ShadowedRicianRandGen(b,m,Omega,N,a)
% This function generates random number according to shadowed Rician
% density function.
%
% INPUTS:
%           b = Scalar (real), Average power of multipath component
%           m = Scalar (real), Fading severity parameter
%       Omega = Scalar (real), Average power of LOS component
%           N = Scalar (real) specifying number of random number to be 
%               generated
% OUTPUTS:
%           X = Scalar (Column Vector if N > 1) specifying random number
%               generated using Shadowed Rician distribution function
% 
% USAGE EXAMPLES:
% X = ShadowedRicianRandGen(0.279,2,0.251);
% 
% REFERENCES:
% A. Abdi, W. C. Lau, M.-S. Alouini, and M. Kaveh, “A new simple model
% for land mobile satellite channels: First- and second-order statistics,”
% IEEE Trans. Wireless Commun., vol. 2, no. 3, pp. 519–528, May 2003.
% Jeruchim, M. C., P. Balaban, and K. S. Shanmugam, Simulation of 
% Communication Systems, New York, Plenum Press, 1992.
%
% Implemented By:
% Ashish (MEET) Meshram
% meetashish85@gmail.com;


% Checking Input Arguments
if nargin<5||isempty(a),a = 10;end
if nargin<4||isempty(N),N = 10000;end
if nargin<3||isempty(Omega)
    error('Missing Input Argument: Please specify omega');
end
if nargin<2||isempty(m)
    error('Missing Input Argument: Please specify m');
end
if nargin<1||isempty(b)
    error('Missing Input Argument: Please specify b');
end


% Implementation Starts Here
X = zeros(N,1);                   % Preallocating memory space for X

% Intermediate Variables 
alpha = ((2*b*m)/(2*b*m + Omega))^m;
beta = Omega/(2*b*(2*b*m + Omega));
lambda = 1/(2*b);


% Maximum value of Shadowed Rician value occurs at x = 0;
maxfx = alpha*lambda;
c = maxfx;
% Accept and Reject Algorithm
for k = 1:N
    accept = false;
    while accept == false
        U2 = c*rand;              % Generating U2, Uniformly disributed 
                                  % random number [0,c]
        U1 = a*rand;              % Generating U1, Uniformly distributed
                                  % in [0,a]
        % Evaluating fx for U1                        
        fx = alpha*lambda*exp(-U1*lambda)*Kummer(m,1,beta*U1);
        % if U2 is less than or equal to fx at U1 then its taken as X else
        % repeat the above procedure
        if U2 <= fx
            X(k) = U1;
            accept = true;
        end
    end
end