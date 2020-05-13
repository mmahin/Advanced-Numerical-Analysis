%{
---------------------------------------------------
 This file uses matlab to calculate Dual Support_VectorsM
 for different parameter combination of C and Gamma
 ---------------------------------------------------
 input > data.txt (feature 1, feature 2, class)
 output > Plot Dual SVM with decision Boundary
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
% Define a 1*n size vector where all values are 1
onevec  = ones(1,n);
% Define a 1*n size vector where all values are 0
zerovec = zeros(1,n);
% Retrieve the labels for Input Data
Labels  = (GetData(:,end))';
Y       = diag(Labels);
% Get the value for y (last column)
y       = (GetData(:,end));


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
        % Create a 1*n vector each value containing the value of C
        C_vector = C_i*onevec;
        % Initialize n*n Gram Matrix containing zero values
        G = zeros(n,n);
        % Define beq for Quadratic Solver
        beq = (0);

        % Gram Matrix Calculation using the Kernel
        for i = 1:n
            for j = 1:n
                G(i,j) = exp(-1*Gamma_i*((norm(Inputs(i,:) - Inputs(j,:)))^2));  
            end
        end

        % Calculate YGY
        YGY = Y*G*Y;

        % Solve the quadratic problem
        [alpha,func_value,exitflag,output,lambda] = quadprog(YGY, ...
            (-1)*onevec, [], [], Labels, beq, zerovec, C_vector);

        % Print Objective Value
        fprintf('The Objective Value is = %f\n',func_value)


        % Retrieve features from Inputs
        % Feature 1 is the first column
        Feature_1 = Inputs(:,1);
        % Feature 2 is the second column
        Feature_2 = Inputs(:,2);



        %{
        -------------------------------------------------------
        Answer to the question number 1b(i)
        -------------------------------------------------------
        %}

        % Draw Figures 
        figure
        % Write title for the figures
        title(['For Parameters- C=', num2str(C_i), ' and Gamma=', num2str(Gamma_i)])
        % Draw Grids
        grid on
        hold on
        % Scatter Plot Feature 1 (+1) in `red`
        scatter(Feature_1(Labels == 1), Feature_2(Labels == 1), 'r', '+')
        hold on
        % Scatter Plot Feature 2 (-1) in `green`
        scatter(Feature_1(Labels == -1), Feature_2(Labels == -1), 'g', '*')
        hold on
        % Scatter Plot for errors in `blue`
        scatter(Feature_1(alpha > Max_error), Feature_2(alpha > Max_error), 'b', 'o')
        hold on




        %{
        -------------------------------------------------------
        Answer to the question number 1b(ii)
        -------------------------------------------------------
        %}

        % Calculate ?_0
        %---------------------------------
        target = 0;
        for j = 1: n
            if alpha(j)>0 && alpha(j)<C_i
                target = j;
                break;
            end
        end
        kernel=0;
        for k = 1: n
             kernel = kernel + y(k)*alpha(k)* exp(-1*Gamma_i*((norm(Inputs(k,:) - Inputs(target,:)))^2));
        end
        beta0 = 1/(y(target)) - kernel;
        %----------------------------------

        % plot the contour: 
        % Define the Meshgrid using linespace
        [XX,YY] = meshgrid(linspace(0,1), linspace(0,1));
        % Get the size
        [mm,nn] = size(XX);
        % Initiallize ZZ of same size with zero values
        ZZ = zeros(size(XX));

        % Calculate ZZ using the RBF Kernel
        for i=1:mm
            for j=1:nn

                % Initialize sum = 0 
                sum = 0;
                % Define grid as vector
                vec = [XX(i,j),YY(i,j)];

                % Calculate ?(?_i * y_i * exp(- ?_i * K(x_i,x_j)))
                for k = 1: n
                    sum =sum +y(k)*alpha(k)* exp(-1*Gamma_i*((norm(Inputs(k,:) - vec))^2));
                end

                % Calculate the final ZZ as 
                % ?(?_i * y_i * exp(- ?_i * K(x_i,x_j))) + ?_0
                % Where ?_0 = lambda.eqlin
                ZZ(i,j) = sum + beta0;
            end
        end


        % Draw Figures 
        figure
        % Write title for the figures
        title(['For Parameters- C=', num2str(C_i), ' and Gamma=', num2str(Gamma_i)])
        % Draw Grids
        grid on
        hold on

        % Plot the contour
        contour(XX,YY,ZZ,[-1,0,1], 'ShowText', 'on');
        hold on

        % Scatter Plot Feature 1 (+1) in `red`
        scatter(Feature_1(Labels == 1), Feature_2(Labels == 1), 'r', '+')
        hold on
        % Scatter Plot Feature 2 (-1) in `green`
        scatter(Feature_1(Labels == -1), Feature_2(Labels == -1), 'g', '*')
        hold on
        % Scatter Plot for errors in `blue`
        scatter(Feature_1(alpha > Max_error), Feature_2(alpha > Max_error), 'b', 'o')
        hold on



        %{
        -------------------------------------------------------
        Answer to the question number 1b(iii)
        -------------------------------------------------------
        %}

        % Calculate Number of support vectors and Bias
        % Define Beta as vector of size 1*p that contains 0
        Beta = zeros(1,p);

        % Calculate Beta values
        for i = 1:n
           Beta = Beta + alpha(i)*Labels(i)*Inputs(i,:); 
        end

        % Find number of Support Vectors
        % Consider Support Vectors for `alpha > Max_error`
        Support_Vectors = find(alpha > Max_error);
        % Calculate the number of Support Vectors from Size
        [total_Support_Vectors,c] = size(Support_Vectors);

        % Print Total Support vectors
        fprintf('For Parameters- C= %.2f and Gamma= %d\n', C_i, Gamma_i)
        fprintf('Number of Support Vectors = %d\n', total_Support_Vectors)

        % Calculate bias for solutions
        bias = Labels(Support_Vectors(1)) - Beta*Inputs(Support_Vectors(1),:)';
        % Print the Bias
        fprintf('The retrieved bias is %f\n\n', bias)
    end
end


