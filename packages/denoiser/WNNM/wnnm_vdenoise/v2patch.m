function [patchmat,nsigmamat] = v2patch(estv,nosv,opt)
%V2PATCH Splite the whole video/volume into overlapped patches.
%   patchmat=V2PATCH(im,opt) returns the overlapped patches patchmat,
%   where v is a video, opt is the options for the function, and patchmat
%   is a two dimensional matrix with each column corresponds to a patch.
%   See also WNNM_VDENOISE, PATCH2V.
p = opt.patchsize;
[nrow,ncol,nframe] = size(estv);
% total number of patches
totpatchnum = (nrow-p+1)*(ncol-p+1);

patchmat = zeros([p*p totpatchnum nframe],'single');
nospatchmat = zeros([p*p totpatchnum nframe],'single'); % noisy patch matrix (estimate noise level)

% the extraction of patches is performed block-wise, not patch-wise and is
% more efficient.
for iframe = 1:nframe
    b = 1;
    for i = 1:p
        for j = 1:p
            block = estv(i:end-p+i,j:end-p+j,iframe); % all the same block as a row
            nosblock = nosv(i:end-p+i,j:end-p+j,iframe); % all the same block as a row
            patchmat(b,:,iframe) = block(:)';
            nospatchmat(b,:,iframe) = nosblock(:)';
            b = b+1;
        end
    end
end
% estimate the noise level of each patch
nsigmamat = opt.vlambda*sqrt(abs(bsxfun(@minus,squeeze(mean((patchmat-nospatchmat).^2,1)),opt.nsigma^2)));

end

