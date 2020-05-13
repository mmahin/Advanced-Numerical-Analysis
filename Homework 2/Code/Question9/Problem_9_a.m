clc
clear all
m = 50;

for i = 1:5
    %Creating random U,Sigma and V
    [U,] = qr(randn(m));
    [V,] = qr(randn(m));
    S=diag(flipud(sort(rand(50,1))));
    %Calculating the matrix
    A = U*S*V';

    %Calculating the SVD of A
    [U2,S2,V2]= svd(A);
    A2=(U2*S2*V2');
    %Calculating the norm of differences
    diffU(i) = norm(U-U2);
    diffV(i) = norm(V-V2);
    diffS(i) = norm(S-S2);
    diffA(i) = norm(A-A2)
    
end

diffU
diffV
diffS
diffA

hold on
xlabel('Values')
ylabel('UError, VError')
plot(diag(U2'*U))
plot(diag(V2'*V))
%plot(diag(S2'*S))
%plot(diag(A2'*A))
legend('U Eror', 'V Error', 'Sigma Error', 'Norm of A')

hold off



