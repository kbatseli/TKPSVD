function [q,Q]=getQ(dim,n)
% [q,Q]=getQ(dim,n)
% ------------------
% Returns the permutation vector (/matrix) q (/Q) such that Q*vec(A)=vec(Atilde) or
% alternatively vec(A)(q)=vec(Atilde).
%
% q		    =   vector, vector of indices that corresponds with the
%				permutation of vec(A) into vec(Atilde),
%
% Q         =   matrix, permutation matrix that corresponds with the
%				permutation of vec(A) into vec(Atilde),
%
% d 		=	vector, vector containing dimensions of the original A
%				tensor, viz. d=size(A),
%
% n         =   vector, contains the dimension of each of the Kronecker 
%               product factors, e.g. decomposition of a 3-way tensor
% 				into a Kronecker product of 2x2x2 with a 3x2x1 factors is
%				specified by n=[3 2 1 2 2 2].
%
% Reference
% ---------
%
% A constructive arbitrary-degree Kronecker product decomposition of
% tensors
% http://arxiv.org/abs/1317261
%
% 2015-2016, Kim Batselier, Ngai Wong

A=reshape([1:prod(dim)],dim);
d=length(dim);
numfac=length(n)/d;

% first reshape A into a length(n)-way tensor
noe=[];
for i=1:d
    noe=[noe n(i:d:length(n))];
end
A=reshape(A,noe);

% permute A
permutedI=zeros(1,length(n));
for i=1:numfac
    permutedI((i-1)*d+1:i*d)=i:numfac:d*numfac;
end
A=permute(A,permutedI);

% reshape A back into a numfac-way tensor
n=noe(permutedI);
n2=zeros(1,length(n)/d);
for i=1:length(n)/d
    n2(i)=prod(n((i-1)*d+1:i*d));
end
q=A(:);
I=eye(prod(dim));
Q=I(q,:);

end
