function [v,f] = patch2v(patmat,frqmat,dim,patchsize)
%PATCH2V Aggregate overlapped patches to a whole video.
%   [v,f]=PATCH2V(patchmat,frqmat,dim,patchsize) returns the aggregated video
%   from overlapped patches patchmat and the corresponding freqency matrix
%   for each patch frqmat. This function is an inverse of the im2patch
%   function and would preform aggregation block-wise as well.
%   See also WNNM_VDENOISE, V2PATCH. 
p = patchsize;
rownum = dim(1)-p+1;
colnum = dim(2)-p+1;
nframe = size(patmat,3);
rows = [1:rownum];
cols = [1:colnum];

v = zeros(dim); 
f = zeros(dim); % frequency of each pixel in the image appears in the patches

% the aggregation of patches is performed block-wise, not patch-wise and is
% more efficient.
for iframe = 1:nframe
    b = 1;
    cpatmat = patmat(:,:,iframe);
    cfrqmat = frqmat(:,:,iframe);
    for i = 1:p
        for j = 1:p
            v(rows-1+i,cols-1+j,iframe) = v(rows-1+i,cols-1+j,iframe)+...
                reshape(cpatmat(b,:)',[rownum colnum]);
            f(rows-1+i,cols-1+j,iframe) = f(rows-1+i,cols-1+j,iframe)+...
                reshape(cfrqmat(b,:)',[rownum colnum]);
            b = b+1;
        end
    end
end
% v = v./(f+eps);

end

