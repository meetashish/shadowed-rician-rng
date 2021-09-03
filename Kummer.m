function y = Kummer(a,b,z,maxit)
% This function implements 1F1(.;.;.), Confluent Hypergeometric function.
%
% INPUTS:
%       a = Scalar and complex
%       b = Scalar and complex
%       z = Scalar and complex
%   maxit = Scalar and real number specifying maximum number of iteration.
%           Default, maxit = 5;
%
% OUTPUT:
%       y = Scalar and complex
%
% Implemented By:
% Ashish (MEET) Meshram
% meetashish85@gmail.com;

% Checking Input Arguments
if nargin<1||isempty(a)
    error('Missing Input Argument: Please specify a');
end
if nargin<2||isempty(b)
    error('Missing Input Argument: Please specify b');
end
if nargin<3||isempty(z)
    error('Missing Input Argument: Please specify z');
end
if nargin<4||isempty(maxit),maxit = 5;end

% Implementation
ytemp = 1;
for k = 1:maxit
    ytemp = ytemp...
            + PochhammerSymbol(a,k)/(PochhammerSymbol(b,k)...
            * factorial(k))*z^k;
    y = ytemp;
end


function y = PochhammerSymbol(x,n)

if n == 0
    y = 1;
else
    y = 1;
    for k = 1:n
        y = y*(x + k - 1);
    end
end