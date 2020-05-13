clc
clear all
%read file
data = csvread('toy.csv');

[m,n] = size(data)
%convert the data into one dimention
new_data = reshape(data,[m*n,1]);

cvx_begin
    lambda = 1;
    variable A(m*n,1);
    A2 = reshape(A,[m,n]);
    %Calculate difference between horizontal pixels
    theta_row = A2(1:m-1,:) - A2(2:m,:);
   %convert to 1D and calculate 2norm
    theta_row_1D_norm = norm(reshape(theta_row,[(m-1)*n,1]),2);
    %Calculate difference between vertical pixels
    theta_col = A2(:,1:n-1) - A2(:,2:n);
    %convert to 1D and calculate 2norm
    theta_col_1D_norm = norm(reshape(theta_col, [m*(n-1),1]),2);
    % minimize the equation. To average vertical and horizontal function,
    % their sum is divided by two
    minimize(0.5*sum_square(A-new_data) + lambda *( (theta_row_1D_norm + theta_col_1D_norm)/2 ));
    
cvx_end

%Calculate the result
result = reshape(A,[m,n]);
figure('NumberTitle', 'off', 'Name', "2 Norm Image Optimal Value " + cvx_optval)
subplot(1,2,1), imshow(data)
title('Original Image.csv')
subplot(1,2,2), imshow(result)
title('New Image.csv')