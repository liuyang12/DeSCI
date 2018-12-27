function im = patch2im(patchmat,frqmat,dim,patchsize)
%PATCH2IM Aggregate overlapped patches to a whole image.
%   im=PATCH2IM(patchmat,frqmat,dim,patchsize) returns the aggregated image
%   from overlapped patches patchmat and the corresponding freqency matrix
%   for each patch frqmat. This function is an inverse of the im2patch
%   function and would preform aggregation block-wise as well.
%   See also WNNM_IMDENOISE, IM2PATCH. 
p = patchsize;
rownum = dim(1)-p+1;
colnum = dim(2)-p+1;
rows = [1:rownum];
cols = [1:colnum];

im = zeros(dim); 
frq = zeros(dim); % frequency of each pixel in the image appears in the patches

% the aggregation of patches is performed block-wise, not patch-wise and is
% more efficient.
b = 1;
for i = 1:p
    for j = 1:p
        im(rows-1+i,cols-1+j) = im(rows-1+i,cols-1+j)+reshape(patchmat(b,:)',[rownum colnum]);
        frq(rows-1+i,cols-1+j) = frq(rows-1+i,cols-1+j)+reshape(frqmat(b,:)',[rownum colnum]);
        b = b+1;
    end
end
im = im./(frq+eps);

end

