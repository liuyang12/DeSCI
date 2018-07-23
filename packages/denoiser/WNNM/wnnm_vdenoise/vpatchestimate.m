function [estpatchmat,frqpatchmat] = vpatchestimate(cframe,nonlocalarr,vindarr,rawpatchmat,nsigmamat,selfindarr,para)
%VPATCHESTIMATE Estimation of the patches in each non-local patch group.
%   [estpatchmat,frqpatchmat]=VPATCHESTIMATE(nonlocalarr,rawpatchmat,
%   selfindarr,para) returns the estimated patch matrix estpatchmat and
%   the frequency matrix of each pixel in a patch.
%   See also WNNM_VDENOISE, WNNM.

% nframe = size(rawpatchmat,3);
pnum = size(rawpatchmat,2);
% memory allocation of the matrices
estpatchmat = zeros(size(rawpatchmat)); % estimated patch matrix (column-by-column)
frqpatchmat = zeros(size(rawpatchmat)); % frequency matrix of each pixel in a patch (column-by-column)

% estimate each patch successively
% todo: adapt the function to parfor for improved performance
% backup of version without parfor
for i = 1:length(selfindarr)
    rawgroup = rawpatchmat(:,nonlocalarr(:,i)+uint32(vindarr(:,i)-1)*pnum); % raw patch group
    mean_rawgroup = repmat(mean(rawgroup,2),[1 para.patchnum]); % mean patch of the patch group
    rawgroup = rawgroup-mean_rawgroup; % deviation of the raw patch group
    
    % WNNM for each patch group by exploiting low-rank property
    if isfield(para,'nnm') && para.nnm % force NNM for reconstruction
        estgroup = nnm(rawgroup,1);
    else
        estgroup = wnnm(rawgroup,para.c,nsigmamat(selfindarr(i),cframe),para.reweighit);
    end
    estgroup = estgroup+mean_rawgroup;
    
%     % Robust PCA estimation
%     estgroup = rpca(rawgroup,ones(size(rawgroup)),nsigmamat(selfindarr(i)));
    
    estpatchmat(:,nonlocalarr(:,i)+uint32(vindarr(:,i)-1)*pnum) ...
        = estpatchmat(:,nonlocalarr(:,i)+uint32(vindarr(:,i)-1)*pnum) + estgroup;
    frqpatchmat(:,nonlocalarr(:,i)+uint32(vindarr(:,i)-1)*pnum) ...
        = frqpatchmat(:,nonlocalarr(:,i)+uint32(vindarr(:,i)-1)*pnum) + 1;
end

% incomplete version with parfor
% patchnum = para.patchnum;
% reweight = para.reweight;
% c        = para.c;
% parfor i = 1:length(selfindarr)
%     patchind = nonlocalarr(:,i)+uint32(vindarr(:,i)-1)*pnum; % patch index
%     selfind  = selfindarr(i);
%     rawgroup = rawpatchmat(:,patchind); % raw patch group
%     mean_rawgroup = repmat(mean(rawgroup,2),[1 patchnum]); % mean patch of the patch group
%     rawgroup = rawgroup-mean_rawgroup; % deviation of the raw patch group
%     
%     % WNNM for each patch group by exploiting low-rank property
%     estgroup = wnnm(rawgroup,c,nsigmamat(selfind,cframe),reweight);
%     estgroup = estgroup+mean_rawgroup;
%     
% %     % Robust PCA estimation
% %     estgroup = rpca(rawgroup,ones(size(rawgroup)),nsigmamat(selfindarr(i)));
%     
%     estpatchmat(:,patchind) = estpatchmat(:,patchind) + estgroup;
%     frqpatchmat(:,patchind) = frqpatchmat(:,patchind) + 1;
% end


end

