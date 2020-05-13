clc
clear all
%read file
X = csvread("X_train.csv");
y = csvread("Y_train.csv");

[n,d] = size(X)

% 6_2C
% Giving weights to the norm, there are eight groups for eighteen features
% and b0 is another group
%parameters are selected based on question
f = zeros(10000,1);
x = ones(19,1);
XX = [ones(n,1) X ];
J = 9;
t = 0.005;
lambda = 0.02;
p = [ 1 1 1 5 6 2 1 1 1];
for k = 1:10000
    %Calculating the gradient
    x = prox_glasso(lambda,t,p, x-t*(1/n*(XX'*((XX*x) - y))));
    %Calculating the function
    g = 1/(2*n)*norm(XX*x-y)^2;
    beta = 1;
    %Calculating function group by group
    for i=1:9
        % We are selecting 
        pi = p(i);
        % All the features belongs to same groups are calculated and the
        % normalizing 
        g = g + lambda*sqrt(pi)*norm(x(beta:beta+pi-1));
        beta = beta + pi;
    end
    f(k) = g;
end
 f_optimal = 49.9649;
 
 %Nop group is considered. 18 feature + beta 0 total 19
 p = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
 f2 = zeros(10000,1);
for k = 1:10000
    %Calculating the gradient
    x = prox_glasso(lambda,t,p, x-t*(1/n*(XX'*((XX*x) - y))));
    %Calculating the function
    g = 1/(2*n)*norm(XX*x-y)^2;
    beta = 1;
    %Calculating function group by group
    for i=1:19
        % We are selecting 
        pi = p(i);
        % All the features belongs to same groups are calculated and the
        % normalizing 
        g = g + lambda*sqrt(pi)*norm(x(beta:beta+pi-1));
        beta = beta + pi;
    end
    f2(k) = g;
end
 f_optimal = 49.9649; 

semilogy(1:10000, (f-f_optimal)/f_optimal, 1:10000,(f2-f_optimal)/f_optimal);
legend('Grouped', 'Without Grouped');