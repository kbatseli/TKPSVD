function S = perfectshuffle(d,n)
% S = perfectshuffle(d,n)
% -------------------------
% Constructs the n^d-by-n^d perfect shuffle matrix S that characterizes a d-way 
% symmetric tensor A of dimenions n.
%
% S         =   perfect shuffle matrix,
%
% d         =   scalar, order of the symmetric tensor A,
%
% n         =   scalar, dimension of each of the modes of A.
%
% Reference
% ---------
%
% A constructive arbitrary-degree Kronecker product decomposition of tensors
% http://arxiv.org/abs/1317261
%
% 2015, Kim Batselier, Ngai Wong

I=eye(n^d);
S=zeros(n^d,n^d);
for i=1:n^(d-1)
   S((i-1)*n+1:i*n,:)=I(i:n^(d-1):n^d,:);
end
