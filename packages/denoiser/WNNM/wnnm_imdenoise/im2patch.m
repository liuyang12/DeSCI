function [patchmat,nsigmamat] = im2patch(estim,nosim,opt)
%IM2PATCH Splite the whole image into overlapped patches.
%   patchmat=IM2PATCH(im,opt) returns the overlapped patches patchmat,
%   where im is an image, opt is the options for the function, and patchmat
%   is a two dimensional matrix with each column corresponds to a patch.
%   See also WNNM_IMDENOISE, PATCH2IM.
p = opt.patchsize;
% total number of patches
totpatchnum = (size(estim,1)-p+1)*(size(estim,2)-p+1);

patchmat = zeros([p*p,totpatchnum],'single');
nospatchmat = zeros([p*p,totpatchnum],'single'); % noisy patch matrix (estimate noise level)

% the extraction of patches is performed block-wise, not patch-wise and is
% more efficient.
b = 1;
for i = 1:p
    for j = 1:p
        block = estim(i:end-p+i,j:end-p+j); % all the same block as a row
        nosblock = nosim(i:end-p+i,j:end-p+j); % all the same block as a row
        patchmat(b,:) = block(:)';
        nospatchmat(b,:) = nosblock(:)';
        b = b+1;
    end
end
% estimate the noise level of each patch
nsigmamat = opt.lambda*sqrt(abs(bsxfun(@minus,mean((patchmat-nospatchmat).^2),opt.nsigma^2)));

end

