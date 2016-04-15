function [B,sigmas]=tkpsvd(A,n,varargin)
% [B,sigmas]=tkpsvd(A,n) or [B,sigmas]=tkpsvd(A,n,R)
% --------------------------------------------------
% Tensor-based Kronecker Product Singular Value Decomposition. Decomposes
% an arbitrary real k-way tensor A into a linear combination of Kronecker
% products as A= \sum_{j=1}^R \sigmas_j A^{dj} \otimes ... \otimes A^{1j},
% with all factors A^{ij} k-way tensors. 
%
% B         =   cell, B{i,j} contains the A^{ij} factor in the TKPSVD,
%
% sigmas    =   vector, contains the coefficients in the linear combination
%               of Kronecker products,
%
% A         =   tensor,
%
% n         =   vector, contains the dimension of each of the Kronecker 
%               product factors, e.g. decomposition of a 3-way tensor
% 				into a Kronecker product of 2x2x2 with a 3x2x1 factors is
%				specified by n=[3 2 1 2 2 2].
%
% R         =   scalar, desired number of terms when calling the sparse
%               version of TTR1SVD.
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
A=reshape(A,n2);

if ~isempty(varargin)
    % R provided, call sparse ttr1
    [U,S,V,sigmas]=ttr1svds(A,varargin{1});
else
    [U,S,V,sigmas]=ttr1svdz(A);
end

% heuristic tolerance to remove numerically zero sigmas
tol=prod(n2)*eps(sigmas(1));
I=find(sigmas>tol);
sigmas(sigmas<tol)=[];
if ~isempty(varargin)
    % R provided, call sparse UV2Uhat
    Uhat=UV2Uhat(U,V,I,varargin{1}*ones(1,length(n2)),n2);
else
     Uhat=UV2Uhat(U,V,I,n2);
end
for j=1:numfac
    for i=1:size(Uhat{1},2)        
        B{j,i}=reshape(Uhat{j}(:,i),n((j-1)*d+1:j*d));
    end
end

end
