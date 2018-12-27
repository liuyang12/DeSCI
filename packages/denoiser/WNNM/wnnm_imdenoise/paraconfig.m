function para = paraconfig(nsigma)
%PARACONFIG Parameter configuration of the WNNM-based image denoiser
%according to the noise level of the noisy image.
%   para=PARACONFIG(nsigma) returns the parameters as a structure for given
%   noise deviation nsigma. Note that those parameters are all empirical 
%   and any adaption of this code is supposed to choose the best ones
%   according to specific applications.
%   See also WNNM_IMDENOISE.
para.nsigma     = nsigma;  % standard deviation of the noisy image
para.windowsize = 30;      % size of the search window
para.delta      = 0.1;     % coefficient of correction between adjacent iterations
  para.adaptboost = true;  % enable adaptive boosting (AB) for WNNM
  para.abeta      = 1;     % eta for noise estimation in WNNM-AB
para.c          = sqrt(2); % constant of weight vector
para.innerloop  = 2;       % numer of inner loops between re-blockmatching
para.reweighit  = 3;       % numer of iterations between re-weight

% choose parameters according to the noise deviation, with larger noise
% deviation accounting for larger patch size, patch number, iteration
% number, and noise estimation constant.
if nsigma<=20     % nsigma<=20
    para.patchsize = 6;    % size of the patch
    para.patchnum  = 70;   % number of patches in a group
    para.iternum   = 8;    % number of iterations
    para.lambda    = 0.54; % constant of noise estimation of each patch
elseif nsigma<=40 % 20<nsigma<=40
    para.patchsize = 7;
    para.patchnum  = 90;
    para.iternum   = 12;
    para.lambda    = 0.56;
elseif nsigma<=60 % 40<nsigma<=60
    para.patchsize = 8;
    para.patchnum  = 120;
    para.iternum   = 14;
    para.lambda    = 0.58;
else              % 60<nsigma
    para.patchsize = 9;
    para.patchnum  = 140;
    para.iternum   = 14;
    para.lambda    = 0.58;
end
% patch step for reduced number of key patches, thus reduced complexity 
% and smaller patch step would hopely result in better performance.
para.patchstep = floor(para.patchsize/2-1); % step of the key patch grid

end
