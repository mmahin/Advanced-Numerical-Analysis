clc
clear all
m = 50;
%format long
for i = 1:5
    %Creating random U,Sigma and V
    [U,] = qr(randn(m));
    [V,] = qr(randn(m));
    S=diag(flipud(sort(rand(50,1))));
    %Calculating the matrix
    A = U*S*V';

    %Calculating the SVD of A
    [U2,S2,V2]= svd(A);
   
    %changing sign column by column
    for j=1:m
        if U2(:,j)'*U(:,j) < 0 
            U2(:,j) = -U2(:,j); 
            V2(:,j) = -V2(:,j);
        end
    end
    %Calculating the norm of differences
    diffU(i) = norm(U-U2);
    diffV(i) = norm(V-V2);
    diffS(i) = norm(S-S2);
    A2=(U2*S2*V2');
    diffA(i) = norm(A-A2);
    con1(i) = cond(A);
    con2(i) = cond(A2);
    
end

diffU
diffV
diffS
diffA
con1
con2


hold on
xlabel('Values')
ylabel('UError, VError')
plot(diag(U2'*U))
plot(diag(V2'*V))
%plot(diag(S2'*S))
%plot(diag(A2'*A))
legend('U Eror', 'V Error')

hold off




