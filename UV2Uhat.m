function Uhat=UV2Uhat(U,V,I,n,varargin)
% Uhat=UV2Uhat(U,V,I,n) or Uhat=UV2Uhat(U,V,I,n,n2) 
% -------------------------------------------------
% Converts the output of the TTR1SVD algorithm into the typical CANDECOMP 
% output as a collection of mode-matrices Uhat.
%
% Uhat       =  cell, Uhat{i} is a matrix containing the mode-i vectors in
%               the CANDECOMP,
%
% U,V 		=  output of the TTr1SVD algorithm,
%
% I 		=	vector, contains indices of nonzero sigmas,
%
% n         =   vector, dimensions of the original tensor or of a tensor
%               that corresponds with a reduced TTr1 tree,
%
% n2        =   vector, dimensions of the original tensor in case n is the
%               dimension of the tensor corresponding with a reduced TTr1
%               tree.
%
% Reference
% ---------
%
% 08/2014, Kim Batselier

indices=leave2ind(I,n);
d=length(n);
Uhat=cell(1,d);
for i=1:d
    if ~isempty(varargin)
        Uhat{i}=zeros(varargin{1}(i),n(i));
    else
        Uhat{i}=zeros(n(i),length(I));
    end
end

for i=1:size(indices,1)
    for j=1:size(indices,2)/2
        Uhat{j}(:,i)=U{indices(i,end-(j-1)*2-1)}(:,indices(i,end-(j-1)*2));
    end
    Uhat{d}(:,i)=V{indices(i,1)}(:,indices(i,2));
end

end
