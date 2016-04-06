%% small demo to illustrate the use of tkpsvd.m

%% centrosymmetric 4-way tensor
clear all;
a=[1:(9*16*4*4)/2 (9*16*4*4)/2:-1:1];
A=reshape(a,[9 16 4 4]);
[Aij,s]=tkpsvd(A,[3 4 2 2 3 4 2 2]);
% all sigmas are distinct, therefore all factors are (skew)-centrosymmetric
s
% inspect centrosymmetry of all factors
norm(fliplr(Aij{1,1}(:))-Aij{1,1}(:))
norm(fliplr(Aij{2,1}(:))-Aij{2,1}(:))
norm(fliplr(Aij{1,2}(:))-Aij{1,2}(:))
norm(fliplr(Aij{2,2}(:))-Aij{2,2}(:))


%% hankel square matrix degree-3 example 
clear all;
clc;
A=hankel([1:27],[27:27+26]);
n=[3 3 3 3 3 3]; 
[B,sigmas]=tkpsvd(A,n);
% inspect Hankel structure of 3rd term
B{3,3}
B{2,3}
B{1,3}

%% symmetric matrix degree 4 example
clear all;
clc;
n=4;
d=4;
A=randn(n^d,n^d);
A=A+A';
dim=n*ones(1,2*d); 
[B,sigmas]=tkpsvd(A,dim);
% inspect (skew)-symmetry structure of 4th term
B{4,4}
B{3,4}
B{2,4}
B{1,4}

