clc
close all
Z = [1,2,3;4,5,6;7,8,7;4,2,3;4,2,4];
%modified gram schidt
[Q1,R1]=mgs(Z)
%householder
[W,R2] = house(Z)
Q2 = formQ(W)
%normal QR
[Q3, R3] = qr(Z,0)
% Difference with main matrix
diff_mgs = norm(Z-(Q1*R1))
diff_house = norm(Z-(Q2*R2))
diff_matlab = norm(Z-(Q3*R3))


%This code is implemented using the sudo code given in the slide. This will
%generate an orthonormal vector Q and triangular vector R from the matrix A
function [Q,R]=mgs(A)
    %get the size of row and collumn in k and l
    [k,l]=size(A);
    %Initialize matrix V having vecotrs v
    V = A;
    %Initialize matrix Q and R
    Q = zeros(k,l);
    R = zeros(l);
    %Calculate Q and R
    for i=1:l
        %Calculate values for the diagonal of the triangular matrix
        R(i,i)= norm(V(:,i));
        %normalize the v vector
        Q(:,i)=V(:,i)/  R(i,i);
        %Calculate other values for R_i vector and modify v to calculate next 
        for j=i+1:l
            R(i,j) = Q(:,i)'*V(:,j);
            V(:,j) = V(:,j)- R(i,j)*Q(:,i);
        end 
    end
    
end


% Here the code is implemented similar to the sudo code written in the
% slide. W will contain v vectors of householder and R will be the
% triangular matrix
function [W,R] = house(A)
%get the size of row and collumn in k and l
[k,l] = size(A);
%Initialize the w to size of A
W = zeros(k,l);
%Calculate the orthonormal and triangular matrix
for i = 1:l
 x = A(i:k,i);
 v = x;
 %calculate the v, only the first value of the vector will be changed. All
 %other will be same
 v(1) = sign(x(1))*norm(x) + x(1);
 v = v/norm(v);
 %A will contain the R matrix
 A(i:k,i:l) = A(i:k,i:l) - 2*v*(v'*A(i:k,i:l));
 %Add v to the w
 W(i:k,i) = v;
end
R = A;
end

% This function will calculate Q frm the W householder vectors. Formula is
% given in the slide
function Q = formQ(W)
%get the size of row and collumn in k and l
[k,l] = size(W);
%Initialize Q with an identity matrix. With this we will reduce the
%opertion of multiplying with I
Q = eye(k);
% Here we will apply householder to Q to calculate Q. We will only consider
% vectors value portion. We will start from the final vector to the first
for i = l:-1:1
 v = W(i:k,i);
 Q(i:k,:) = Q(i:k,:) - 2*v*(v'*Q(i:k,:));
end 
end




