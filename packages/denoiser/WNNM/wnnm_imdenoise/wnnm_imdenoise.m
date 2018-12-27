function [denoisedim,psnrall] = wnnm_imdenoise(noisyim,orgim,para)
%WNNM_IMDENOISE Weighted nuclear norm minimization (WNNM)-based image
%denoising.
%   denoisedim=WNNM_IMDENOISE(noisyim,para) returns the denoised image
%   denoisedim using WNNM, where noisyim is the input noisy image and para
%   is the parameters for the WNNM image denoiser.
%   # Model
%     The basic idea of WNNM for image denoising is that the patch group
%     with nonlocal self-similarity can be expressed as a low-rank matrix.
%     Low-rank matrix approximation problem has already been solved by
%     low-rank matrix factorization (LRMF) and nuclear norm minimization 
%     (NNM). WNNM approach to this problem focuses on varied weights of
%     each singular value, where greater weight indicates less importance.
%     The CVPR'14 paper gives an convergence proof to the fixed point of 
%     the WNNM problem with non-descending weights, which is common in 
%     image denoising applications where small singular values of the patch
%     group are supposed to be surpressed.
%   # Parameters
%     parameters              [default] value
%     para.patchsize
%   # Reference
%   [1]  S. Gu, L. Zhang, W. Zuo, and X. Feng, "Weighted Nuclear Norm 
%          Minimization with Application to Image Denoising," in 2014 IEEE 
%          Conference on Computer Vision and Pattern Recognition (CVPR), 
%          2014, pp. 2862-2869.
%   [2]  S. Gu, Q. Xie, D. Meng, W. Zuo, X. Feng, and L. Zhang, "Weighted 
%          Nuclear Norm Minimization and Its Applications to Low Level 
%          Vision," International Journal of Computer Vision, vol. 121, 
%          no. 2, pp. 183-208, 2017.
%   See also NEIGHBORIND, IM2PATCH, BLOCKMATCH, PATCHESTIMATE, WNNM,
%            PATCH2IM, PARACONFIG.

%   Code copyright(C) 2017-2018 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>. This code is a
%   derivative of WNNM code <a href="matlab: 
%   web('http://www4.comp.polyu.edu.hk/~cslzhang/code/WNNM_code.zip')"> 
%   provided by S. Gu.
%   Last modified Oct 29, 2017.

% [height,width] = size(noisyim);
% patchdim = para.patchsize*para.patchsize; % dimension of the patch (vector)
% totpatchnum ...  % total number of patches in the noisy image
%   = (height-para.patchsize+1)*(width-para.patchsize+1);

% pre-calculation of the indexes of the neighbors within the search window
[neighborindarr,neighbornumarr,selfindarr] = neighborind(size(noisyim),para);

% % [unnecessary] memory allocation of the matrices, which indicates the size
% % of the matrices
% nonlocalarr = zeros(para.patchnum,length(neighbornumarr)); % index array of patches with non-local similarity as the key patch
% rawpatchmat = zeros(patchdim,totpatchnum); % raw patch mattrix (column-by-column)
% estpatchmat = zeros(patchdim,totpatchnum); % estimated patch matrix (column-by-column)
% frqpatchmat = zeros(patchdim,totpatchnum); % frequency matrix of each pixel of a patch (column-by-column)

% WNNM denoisng for several iterations
estim = noisyim;
for iter = 1:para.iternum
    if ~para.adaptboost
        estim = estim + para.delta*(noisyim-estim); % correction between adjacent iterations
    else % adative boosting for WNNM-based image denoising
        Eest = mean(abs(estim(:)).^2);
        En = para.abeta*abs(para.nsigma^2-var(noisyim(:)-estim(:)));
        rho = sqrt(Eest)/(sqrt(Eest)+sqrt(max(Eest-En,0)));
        fprintf('    Iteration % 2d, rho = %.3f. \n',iter,rho);
        estim = estim + (1-rho)*(noisyim-estim); % correction between adjacent iterations
    end
    [rawpatchmat,nsigmamat] = im2patch(estim,noisyim,para); % splite the whole image into overlapped patches
    if mod(iter-1,para.innerloop)==0
        para.patchnum = para.patchnum-10; % less non-local patches with loweer noise level
        % calculate the patches with non-local similarity for each key patch
        nonlocalarr = blockmatch(rawpatchmat,neighborindarr,neighbornumarr,selfindarr,para);
        if iter==1 % initial noise level of each patch
            nsigmamat = para.nsigma*ones(size(nsigmamat));
        end
    end
    % patch estimation by means of WNNM
    [estpatchmat,frqpatchmat] = patchestimate(nonlocalarr,rawpatchmat,nsigmamat,selfindarr,para);
    % aggregate overlapped patches to the whole image
    estim = patch2im(estpatchmat,frqpatchmat,size(noisyim),para.patchsize);
    
    if ~isempty(orgim)
        % error statistics and reconstruction quality accessment
        rmseall(iter) = sqrt(immse(estim,orgim));
        psnrall(iter) = psnr(estim,orgim,255);
        ssimall(iter) = ssim(estim,orgim);
        fprintf('  WNNM iteration % 2d, RMSE = %2.2f, PSNR = %2.2f dB, SSIM = %.4f.\n',iter,rmseall(iter),psnrall(iter),ssimall(iter));
    end
end
denoisedim = estim;

end

