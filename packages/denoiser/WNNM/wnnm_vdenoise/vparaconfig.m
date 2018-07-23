function vpara = vparaconfig(nsigma)
%PARACONFIG Parameter configuration of the WNNM-based image denoiser
%according to the noise level of the noisy image.
%   para=PARACONFIG(nsigma) returns the parameters as a structure for given
%   noise deviation nsigma. Note that those parameters are all empirical 
%   and any adaption of this code is supposed to choose the best ones
%   according to specific applications.
%   See also WNNM_IMDENOISE.
vpara.nsigma     = nsigma;  % standard deviation of the noisy image
vpara.windowsize = 30;      % size of the search window (centered)
vpara.sframesize = 4;       % search frame size (centered)
vpara.delta      = 0.1;     % coefficient of correction between adjacent iterations
  vpara.adaptboost = false; % enable adaptice boosting (AB) in WNNM
  vpara.abeta      = 0.25;  % eta in AB noise estimation
vpara.c          = sqrt(2); % constant of weight vector
vpara.innerloop  = 2;       % numer of inner loops between re-blockmatching
vpara.reweighit  = 3;       % numer of iterations between re-weight

% choose parameters according to the noise deviation, with larger noise
% deviation accounting for larger patch size, patch number, iteration
% number, and noise estimation constant.
if nsigma<=20     % nsigma<=20
    vpara.patchsize = 6;    % size of the patch
    vpara.patchnum  = 70;   % number of patches in a group
    vpara.iternum   = 8;    % number of iterations
    vpara.vlambda   = 0.54; % constant of noise estimation of each patch
elseif nsigma<=40 % 20<nsigma<=40
    vpara.patchsize = 7;
    vpara.patchnum  = 90;
    vpara.iternum   = 12;
    vpara.vlambda   = 0.56;
elseif nsigma<=60 % 40<nsigma<=60
    vpara.patchsize = 8;
    vpara.patchnum  = 120;
    vpara.iternum   = 14;
    vpara.vlambda   = 0.58;
else              % 60<nsigma
    vpara.patchsize = 9;
    vpara.patchnum  = 140;
    vpara.iternum   = 14;
    vpara.vlambda   = 0.58;
end
% patch step for reduced number of key patches, thus reduced complexity 
% and smaller patch step would hopely result in better performance.
vpara.patchstep = floor(vpara.patchsize/2-1); % step of the key patch grid

end
