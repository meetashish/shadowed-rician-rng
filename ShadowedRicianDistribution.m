function [f,F] = ShadowedRicianDistribution(b,m,Omega,x)

lambda = 1/(2*b);
alpha = (2*b*m)/(2*b*m + Omega);
beta = Omega/(2*b*(2*b*m + Omega));

% Theoretical PDF
f = zeros(length(x),1);
for k = 1:length(x)
    f(k) = (alpha^m)*lambda*exp(-x(k)*lambda)*Kummer(m,1,beta*x(k));
end

% Theoretical CDF
sumk = zeros(500,1);
F = zeros(length(x),1);
for p = 1:length(x)
    for q = 1:500
        mmk = gamma(m+k);
        mk = gamma(m);
        betabylambdak = (beta/lambda)^k;
        gammak = gammainc(k+1,lambda*x(p));
        sumk(q) = (mmk/(mk*(factorial(k))^2))*betabylambdak*gammak;
    end
    F(p) = alpha*sum(sumk);
end
