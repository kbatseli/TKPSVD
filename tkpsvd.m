function [B,sigmas]=tkpsvd(A,n,varargin)
% [B,sigmas]=tkpsvd(A,n) or [B,sigmas]=tkpsvd(A,n,R)
% --------------------------------------------------
% Tensor-based Kronecker Product Singular Value Decomposition. Decomposes
% an arbitrary real matrix A into a linear combination of Kronecker
% products as A= \sum_{j=1}^R \sigmas_j A^{dj} \otimes ... \otimes A^{1j}. 
%
% B         =   cell, B{i,j} contains the A^{ij} factor in the TKPSVD,
%
% sigmas    =   vector, contains the coefficients in the linear combination
%               of Kronecker products,
%
% A         =   matrix,
%
% n         =   vector, contains the dimension of each of the Kronecker 
%               as [m_1 n_1 m_2 n_2 .... m_d n_d] with A^{ij} an m_i-by-n_i
%               matrix,
%
% R         =   scalar, desired number of terms when calling the sparse
%               version of TTR1SVD.
%
% Reference
% ---------
%
% A constructive arbitrary-degree Kronecker product decomposition of matrices
%
% 2015, Kim Batselier, Ngai Wong

% n contains the dimensions of the even order tensor n = [n1 n2 n3 ... nd]
d=length(n);
if mod(d,2)~=0
    error('Kronecker Product SVD only works for even order tensors')
end

% reshape matrix A into a d-way tensor
io=1:2:d;   % odd indices
ie=2:2:d;   % even indices   
noe=[n(io) n(ie)];
A=reshape(A,noe);

% permute A
permutedI=zeros(1,d);
permutedI(1:2:end)=1:d/2;
permutedI(2:2:end)=d/2+1:d;
A=permute(A,permutedI);

n=noe(permutedI);
n2=zeros(1,d/2);
for i=1:d/2
    n2(i)=prod(n((i-1)*2+1:i*2));
end
% reshape A into a d/2-way tensor
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
for j=1:d/2
    for i=1:size(Uhat{1},2)        
        B{j,i}=reshape(Uhat{j}(:,i),n((j-1)*2+1:j*2));
    end
end

end
