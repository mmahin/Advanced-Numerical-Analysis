%{
---------------------------------------------------
 This file uses matlab to plot Convergence Histories
 for different parameter combination of C and Gamma
 ---------------------------------------------------
 input > data.txt (feature 1, feature 2, class)
 output > Plot Convergence History
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
% Get value for y from the third column
y  = GetData(:,3);
% Get size of y
m = size(y);

% Create a vector of size 1*n of values 1
onevec  = ones(1,n);
% Create a vector of size 1*n of values 2
zerovec = zeros(1,n);
% Get label from the data file
Labels  = (GetData(:,end))';
Y       = diag(Labels);

% Given Values of Gamma and C
Gamma   = [10 50 100 500];    
C       = [0.01 0.1 0.5 1];          


% for loop to access all combinations
for gamma_val = 1:4
    for c_val = 1:4
    
        % Get parameter value combinations from both lists
        % Select value for Gamma
        Gamma_i = Gamma(gamma_val);
        % Select value for C
        C_i = C(c_val);
        % Create a 1*n vector each value containing the value of C
        C_vector = C_i*onevec;
        % Define beq for Quadratic Solver
        beq = (0);
        
        % Initialize n*n Gram Matrix containing zero values
        G = zeros(n,n);
        
        % Gram Matrix Calculation using the Kernel
        for i = 1:n
            for j = 1:n
                G(i,j) = exp(-1*Gamma_i*((norm(Inputs(i,:) - Inputs(j,:)))^2));  
            end
        end

        % Calculate YGY
        YGY = Y*G*Y;

        
        % Solve the quadratic problem
        [alpha1,func_value,exitflag,output,lambda] = quadprog(YGY, ...
            (-1)*onevec, [], [], Labels, beq, zerovec, C_vector);
        
        % Print Objective Value
        fprintf('The Objective Value is = %f',func_value)
        
        % Define beta
        beta = 0.5;
        % Initiate Alpha of size `y` of zero values
        alpha = zeros(size(y));
        % Consider k=50 (depends on the iteration we require)
        k = 50;
        % Initiate error vector of size k of zero values
        errors = zeros(k);
        K = zeros(k);
        % Iterate P_D(x) function to calculate convergence until error
        for i = 1:k
            t = 1;
            gradient_t = Y*G*Y*alpha - 1;
            % Call the P_D(x) function
            PD = P_D((alpha - (t*gradient_t)),y,C_i);
            G_t = (alpha - PD)/t;
            
            % Call function `f` in loop
            while f((alpha - (t * G_t)),Y,G) > (f(alpha,Y,G) ...
                    - (t * gradient_t' * G_t) + ((t/2)*(norm(G_t))^2))
                t = beta * t;
            end
            % Calculate Alpha
            alpha = P_D((alpha - (t*gradient_t)),y,C_i);
            % Calculate new f^*
            f_star = f(alpha,Y,G) 
            error(i)=(f_star -func_value);
            if (error(i)<Max_error)
                break
            end
        end
        
       
        % Draw Figures 
        figure
        % Write title for the figures
        title(['For Parameters- C=', num2str(C_i), ' and Gamma=', num2str(Gamma_i)])
        % Draw Grids
        grid on
        hold on
        plot(error)
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
    

    
   