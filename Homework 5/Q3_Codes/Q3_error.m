
%{
---------------------------------------------------
 This file uses matlab to calculate normalized output
 difference between CVX and P_D(x) function
 for different parameter combination of C and Gamma
 ---------------------------------------------------
 input > data.txt (feature 1, feature 2, class)
 output > Calculated difference of errors between
            CVX and P_D(x) function
 --------------------------------------------------
 Compile > Open the Source file using Matlab and Run
 --------------------------------------------------
 Author > 
    Md Mahin (PSID: 1900421)
    Shanto Roy (PSID: 1894941)
%}


% Clear All
clc
clear all

% Input data from the given file `data.txt`
% get the data filename
filename = 'data.txt';
% delimiter used in the data file
delimiter = ' ';
% Read data from the file
GetData = importdata(filename, delimiter);

% Maximum tolerable Error
Max_error = 1e-5;

% Get feature 1 and Feature 2 from imported data
Inputs  = GetData(:,1:2);
% Get size of the input data
[n, p]  = size(Inputs);
% Get the value for y (last column)
y  = GetData(:,3);


% lists of given parameters for gamma and C
% List of Gamma
Gamma   = [10 50 100 500];    
% List of C
C       = [0.01 0.1 0.5 1];            


% for loop to access all combinations
for gamma_val = 1:4
    for c_val = 1:4
    
        % Get parameter value combinations from both lists
        % Select value for Gamma
        Gamma_i = Gamma(gamma_val);
        % Select value for C
        C_i = C(c_val);
        
        % Calculate CVX outputs
        for i=1:10
        x = randn(n,1);
        %-------------------CVX Calculation-------------------
        cvx_begin
            variable z(n);
            minimize( 0.5*(norm(z-x))^2 )
            subject to
                z'*y == 0
                0 <= z <= C_i
        cvx_end
        %-----------------------------------------------------
        %--------------------P_D Calculation------------------
        q = P_D(x,y,C_i);
        %-----------------------------------------------------
        % Differentiate CVX and P_D outputs using Normalization
        norm(z-q)
        end
    end
end


% Define the P_D function
    %{
        ---------------Function Starts Here----------------
    %}
    function y = P_D(x,y,C_i)
        Max_error = 1e-5;
        mu1 = 0;
        mu2 = 10000;
        
        %{
        --------------------------------------------------
        Working Code Block, if not used seperate Functions
        The seperate functions are defined below as `f, P_box, Phi`
        --------------------------------------------------
        mu_star = (mu1 + mu2)/2;
        s = Phi(mu_star,x,y,C_i);
        inequality = (mu2-mu1)/2;
        while ((s ~= 0) & (inequality > Max_error))
           s1 = sign(Phi(mu_star,x,y,C_i));
           s2 = sign(Phi(mu1,x,y,C_i));
           if(s1==s2)
               mu1 = mu_star;
           else
               mu2 = mu_star;
           end
           s = Phi(mu_star,x,y,C_i);
           inequality = (mu2-mu1)/2;
           mu_star = (mu1 + mu2)/2;
        end
        mu_star = (mu1 + mu2)/2;
        ---------------------------------------------------
        %}
        
        mu_star = 0;
        N = 1;
        while N < 1000
            mu_star = (mu1 + mu2)/2; 
            if Phi(mu_star,x,y,C_i) == 0 | ((mu2 - mu1)/2) < Max_error 
                break
            end
            N = N + 1; 
             if sign(Phi(mu_star,x,y,C_i)) == sign(Phi(mu1,x,y,C_i))  
                 mu1 = mu_star; 
             else
                 mu2 = mu_star; 
             end
        end 
        % Calls Function `P_box`
        y =P_box((x - (mu_star * y)),C_i);
    
    end
    %-------------------------------------------------------

        
    
    % Define Function `Phi`
    %{
    ---------------Function Starts Here----------------
    %}
    function z = Phi(mu,x,y,C_i)
        % Calls Function `P_box`
        z= y'*P_box((x - (mu * y)),C_i);
    end
    %--------------------------------------------------



    % Define Function P_box
    %{
    ---------------Function Starts Here----------------
    %}
    function w = P_box(q,C_i)
        w = zeros(size(q));
        for i = 1: size(q)
            if q(i) < 0
                w(i) = 0;
            elseif q(i) > C_i
                w(i) = C_i;
            else
                w(i) = q(i);
            end
        end

    end
    %---------------------------------------------------
        
        
        
        
    % Define function `f` that returns f0
    %{
    ---------------Function Starts Here----------------
    %}
    function fo = f(q,Y,G)
        fo = (((1/2)*q'*Y*G*Y*q)-sum(q));
    end
    %--------------------------------------------------

    