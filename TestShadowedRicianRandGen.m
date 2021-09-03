clear;clc;
%% Parameters
N = 1e5;
b = 0.279;
m = 2;
Omega = 0.251;

%% Generate Shadowed Rician Random Number
X = ShadowedRicianRandGen(b,m,Omega,N);

%% Points for which distribution has to be evaluated
x = linspace(0,5,1000);

%% Estimate Distribution
[fsim,Fsim] = EstimateDistribution(X,x);

%% Theoretical PDF & CDF
[fana,Fana] = ShadowedRicianDistribution(b,m,Omega,x);

%% Plot Results
subplot(121);plot(x(1:30:end),Fsim(1:30:end),'-b*',x,Fana,'-');grid on;
xlabel('x');ylabel('F_X(x) = P(X<=x)');
subplot(122);plot(x(1:30:end),fsim(1:30:end),'-b*',x,fana,'-');grid on;
xlabel('x');
ylabel('Probability Density Function');
legend('Simulated','Theoretical');