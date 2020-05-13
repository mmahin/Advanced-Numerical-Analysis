clc
clear all
close all

%Single Precision Number Test

t = 24;
n1 = single(2^t - 3);
n2 = single(2^t - 2);
n3 = single(2^t - 1);
n4 = single(2^t);
n5 = single(2^t + 1);
n6 = single(2^t + 2);
n7 = single(2^t + 3);
%Calculate difference with most adjacent number
diff_32(2) = n3-n2;
diff_32(3) = n4-n3;
diff_32(4) = n5-n4;
diff_32(5) = n6-n5;
diff_32(6) = n7-n6;
%If difference is 1 it exists, if 0/2 it does not exist
diff_32

%Double Precision Number Test
t = 53;

n1 = double(2^t - 3);
n2 = double(2^t - 2);
n3 = double(2^t - 1);
n4 = double(2^t);
n5 = double(2^t + 1);
n6 = double(2^t + 2);
n7 = double(2^t + 3);
%Calculate difference with most adjacent number
diff_64(1) = n2-n1;
diff_64(2) = n3-n2;
diff_64(3) = n4-n3;
diff_64(4) = n5-n4;
diff_64(5) = n6-n5;
diff_64(6) = n7-n6;
%If difference is 1 it exists, if 0/2 it does not exist
diff_64