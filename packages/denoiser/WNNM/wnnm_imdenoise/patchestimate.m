function [estpatchmat,frqpatchmat] = patchestimate(nonlocalarr,rawpatchmat,nsigmamat,selfindarr,para)
%PATCHESTIMATE Estimation of the patches in each non-local patch group.
%   [estpatchmat,frqpatchmat]=PATCHESTIMATE(nonlocalarr,rawpatchmat,
%   selfindarr,para) returns the estimated patch matrix estpatchmat and
%   the frequency matrix of each pixel in a patch.
%   See also WNNM_IMDENOISE, WNNM.

% memory allocation of the matrices
estpatchmat = zeros(size(rawpatchmat)); % estimated patch matrix (column-by-column)
frqpatchmat = zeros(size(rawpatchmat)); % frequency matrix of each pixel in a patch (column-by-column)

% estimate each patch successively
for i = 1:length(selfindarr)
    rawgroup = rawpatchmat(:,nonlocalarr(1:para.patchnum,i)); % raw patch group
    mean_rawgroup = repmat(mean(rawgroup,2),[1 para.patchnum]); % mean patch of the patch group
    rawgroup = rawgroup-mean_rawgroup; % deviation of the raw patch group
    
    % WNNM for each patch group by exploiting low-rank property
    estgroup = wnnm(rawgroup,para.c,nsigmamat(selfindarr(i)),para.reweighit);
    estgroup = estgroup+mean_rawgroup;
    
%     % Robust PCA estimation
%     estgroup = rpca(rawgroup,ones(size(rawgroup)),nsigmamat(selfindarr(i)));
    
    estpatchmat(:,nonlocalarr(1:para.patchnum,i)) ...
        = estpatchmat(:,nonlocalarr(1:para.patchnum,i)) + estgroup;
    frqpatchmat(:,nonlocalarr(1:para.patchnum,i)) ...
        = bsxfun(@plus,1,frqpatchmat(:,nonlocalarr(1:para.patchnum,i)));
end

end

