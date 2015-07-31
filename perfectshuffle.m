function S = perfectshuffle(p,q)
% S = perfectshuffle(p,q)
% -------------------------
% Constructs the pq-by-pq perfect shuffle matrix S for a p-by-q matrix A.
%
% S         =   perfect shuffle matrix,
%
% p         =   scalar, number of rows of A,
%
% q         =   scalar, number of columns of A.
%
% Reference
% ---------
%
% A constructive arbitrary-degree Kronecker product decomposition of matrices
% http://arxiv.org/abs/1317261
%
% 2015, Kim Batselier, Ngai Wong

r=p*q;
I=eye(r);
S=[];
for i=1:q
   S=[S; I(i:q:r,:)];
end
