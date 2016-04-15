function Atilde=getAtilde(A,n)
% Atilde=getAtilde(A,n)
% ---------------------
% Returns the permuted and reshaped tensor A from which an orthogonal
% polyadic decomposition needs to be computed in order to compute its 
% TKPSVD as specified by the vector n.
%
% Atilde    =   tensor, permuted and reshaped A,
%
% A         =   tensor, multi-way array from which the TKPSVD needs to
%				be computed,
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

dim=size(A);
d=length(dim);
numfac=length(n)/d;
if rem(length(n),d) ~=0
    error('Not all factors have the same order.')
end

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
Atilde=reshape(A,n2);
end
