function C = tkron(A,B)
%  C = tkron(A,B)
% ---------------
% Returns the Kronecker product of two input tensors A,B.
%
% C         =   multiway array, Kronecker product of A with B.
%
% A,B       =   multiway array, tensors of the same order.
%
% Reference
% ---------
%
% 01/206, Zhongming Chen

% C{i,j,k} = C{(i2,i1), (j2,j1), (k2, k1)} = B{i2,j2,k2} * A{i1,j1,k1}
siz1 = size(A);   siz2 = size(B);
if length(siz1) > length(siz2)
    siz2=[siz2 ones(1,length(siz1)-length(siz2))];
else
    siz1=[siz1 ones(1,length(siz2)-length(siz1))];
end
d = length(siz1);    siz = siz1.*siz2;
N = prod(siz);       C = zeros(siz);
for ndx = 1:N
    ind = ind2sub(siz, ndx);
    ind1 = zeros(1,d);   ind2 = zeros(1,d);
    for t = 1:d
        temp_siz = [siz2(t) siz1(t)];        % Be careful
        temp_ind = ind2sub(temp_siz, ind(t));
        ind1(t) = temp_ind(2);   ind2(t) = temp_ind(1);  % Be careful
    end
    ndx1 = sub2ind(siz1, ind1);
    ndx2 = sub2ind(siz2, ind2);
    C(ndx) = A(ndx1)*B(ndx2);
end

end


function ind = ind2sub(siz, ndx)
nout = size(siz,2);
ind=zeros(1,nout);
if length(siz)<=nout,
   siz = [siz ones(1,nout-length(siz))];
 else
   siz = [siz(1:nout-1) prod(siz(nout:end))];
 end
 n = length(siz);
 k = [1 cumprod(siz(1:end-1))];
 for i = n:-1:1,
   vi = rem(ndx-1, k(i)) + 1;         
   vj = (ndx - vi)/k(i) + 1; 
   ind(i) = vj; 
   ndx = vi;     
 end
end


function ndx = sub2ind(siz, ind)
a = ind - ones(size(ind));   a(1) = a(1) + 1;
b = [1 cumprod(siz(1:end-1))];
ndx = a * b';

end