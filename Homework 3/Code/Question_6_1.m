clc
clear all
%read file
X = csvread("X_train.csv");
y = csvread("Y_train.csv");
XX = [ones(m,1) X];
[m,n] = size(X)

 cvx_begin
    lambda = 1;
     variable a(n+1)
     %equation
     minimize( sum_square(XX*a-y)/(2*m) + lambda/2*sum_square(a(2:19)) )
 cvx_end



