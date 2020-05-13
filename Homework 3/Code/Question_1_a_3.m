clc
clear all
num = csvread('baboon.csv');
% imshow(num)
% imhist(num,100)
clc
clear all
%read file
data = csvread('baboon.csv');

[m,n] = size(data)
%convert the data into one dimention
new_data = reshape(data,[m*n,1]);

for k=0:8
    cvx_begin
    lambda = 10^(-k/4);
    variable A(m*n,1);
    A2 = reshape(A,[m,n]);
    %Calculate difference between horizontal pixels
    theta_row = A2(1:m-1,:) - A2(2:m,:);
   %convert to 1D and calculate norm
    theta_row_1D_norm = norm(reshape(theta_row,[(m-1)*n,1]),1);
    %Calculate difference between vertical pixels
    theta_col = A2(:,1:n-1) - A2(:,2:n);
    %convert to 1D and calculate norm
    theta_col_1D_norm = norm(reshape(theta_col, [m*(n-1),1]),1);
    % minimize the equation. To average vertical and horizontal function,
    % their sum is divided by two
    minimize(0.5*sum_square(A-new_data) + lambda *( (theta_row_1D_norm + theta_col_1D_norm)/2 ));
    
    cvx_end
    
    
    
    %Calculate the result
    result1 = reshape(A,[m,n]);
    figure('NumberTitle', 'off', 'Name', 'Anisotropic 1 Norm Image' + "k = " + k + " Optimal Value " + cvx_optval)
    subplot(1,2,1), imshow(data)
    title('Original Image.csv')
    subplot(1,2,2), imshow(result1)
    title('New Image.csv')
    
    figure('NumberTitle', 'off', 'Name', 'Anisotropic 1 Norm Image' + "k = " + k + " Optimal Value " + cvx_optval)
    subplot(1,2,1), histogram(data, 100, 'BinLimits', [0,1])
    title('Original Image.csv')
    subplot(1,2,2, 'align'), histogram(result1, 100, 'BinLimits', [0,1])
    title('New Image.csv')
 
        cvx_begin
    lambda = 10^(-k/4);
    variable A(m*n,1);
    A2 = reshape(A,[m,n]);
    %Calculate difference between horizontal pixels
    theta_row = A2(1:m-1,:) - A2(2:m,:);
   %convert to 1D and calculate 2 norm
    theta_row_1D_norm = norm(reshape(theta_row,[(m-1)*n,1]),2);
    %Calculate difference between vertical pixels
    theta_col = A2(:,1:n-1) - A2(:,2:n);
    %convert to 1D and calculate 2 norm
    theta_col_1D_norm = norm(reshape(theta_col, [m*(n-1),1]),2);
    % minimize the equation. To average vertical and horizontal function,
    % their sum is divided by two
    minimize(0.5*sum_square(A-new_data) + lambda *( (theta_row_1D_norm + theta_col_1D_norm)/2 ));
    
    cvx_end
    
    
    
    %Calculate the result
    result2 = reshape(A,[m,n]);
    figure('NumberTitle', 'off', 'Name', 'Isotropic 2 Norm Image' + "k = " + k + " Optimal Value " + cvx_optval)
    subplot(1,2,1), imshow(data)
    title('Original Image.csv')
    subplot(1,2,2), imshow(result2)
    title('New Image.csv')
    
    figure('NumberTitle', 'off', 'Name', 'Isotropic 2 Norm Image' + "k = " + k + " Optimal Value " + cvx_optval)
    subplot(1,2,1), histogram(data, 100, 'BinLimits', [0,1])
    title('Original Image.csv')
    subplot(1,2,2, 'align'), histogram(result2, 100, 'BinLimits', [0,1])
    title('New Image.csv')
end


