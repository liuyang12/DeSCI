%TEST_DESCI_CASSI_BIRD Test decompress snapshot compressive imaging (DeSCI) 
%for simulated coded aperture snapshot spectral imaging (CASSI) `bird` 
%dataset.
% Reference
%   [1] Y. Liu, X. Yuan, J. Suo, D.J. Brady, and Q. Dai, Rank Minimization 
%       for Snapshot Compressive Imaging, IEEE Trans. Pattern Anal. Mach. 
%       Intell. (TPAMI), DOI:10.1109/TPAMI.2018.2873587, 2018.
%   [2] X. Yuan, Generalized alternating projection based total variation 
%       minimization for compressive sensing, in Proc. IEEE Int. Conf. 
%       Image Process. (ICIP), pp. 2539-2543, 2016.
% Dataset
%   `toy` from the CAVE multispectral image database
%     (http://www1.cs.columbia.edu/CAVE/databases/multispectral/).
%   `bird` from the multiframe CASSI system [3] captured in [4].
%     [3] D. Kittle, K. Choi, A. Wagadarikar, and D. J. Brady, Multiframe 
%         image estimation for coded aperture snapshot spectral imagers, 
%         Appl. Opt., vol. 49, no. 36, pp. 6824-6833, 2010.
%     [4] A. Rajwade, D. Kittle, T.-H. Tsai, D. Brady, and L. Carin, Coded 
%         Hyperspectral Imaging and Blind Compressive Sensing, SIAM J. on 
%         Imag. Sci., vol. 6, no. 2, pp. 782-812, 2013.
% Contact
%   Xin Yuan, Bell Labs, xyuan@bell-labs.com, initial version Jul 2, 2015.
%   Yang Liu, Tsinghua University, y-liu16@mails.tsinghua.edu.cn, last 
%     update Dec 17, 2018.
%   See also GAPDENOISE_CACTI, GAPDENOISE.
clear; clc;
% [0] environment configuration
addpath(genpath('../algorithms')); % algorithms
addpath(genpath('../packages')); % packages
addpath(genpath('../utils')); % utilities

datasetdir = '../dataset'; % dataset
resultdir  = '../results'; % results

% [1] load dataset
para.type   = 'cassi'; % type of dataset, cassi or cacti
para.name   = 'bird'; % name of dataset
para.number = 24; % number of frames in the dataset

datapath = sprintf('%s/%s%d_%s.mat',datasetdir,para.name,...
    para.number,para.type);

if exist(datapath,'file')
    load(datapath); % mask, meas, orig (and para)
else
    error('File %s does not exist, please check dataset directory!',...
        datapath);
end

para.nframe = 1; % number of coded frames in this test
para.MAXB   = 255;

[nrow,ncol,nmask] = size(mask);
nframe = para.nframe; % number of coded frames in this test
MAXB = para.MAXB;

% [2] apply GAP-Denoise for reconstruction
% yall = meas/MAXB;

para.Mfunc  = @(z) A_xy(z,mask);
para.Mtfunc = @(z) At_xy_nonorm(z,mask);

para.Phisum = sum(mask.^2,3);
para.Phisum(para.Phisum==0) = 1;
% common parameters
para.lambda   =     1; % correction coefficiency
para.acc      =     1; % enable GAP-acceleration
para.flag_iqa = false; % disable image quality assessments in iterations

%% [2.1] GAP-TV, ICIP'16
para.denoiser = 'tv'; % TV denoising
  para.maxiter  = 250; % maximum iteration
  para.tvweight = 5; % weight for TV denoising
  para.tviter   = 5; % number of iteration for TV denoising
  
[vgaptv,psnr_gaptv,ssim_gaptv,tgaptv] = ...
    gapdenoise_cacti(mask,meas,orig,[],para);

fprintf('GAP-%s mean PSNR %2.2f dB, mean SSIM %.4f, total time % 4.1f s.\n',...
    upper(para.denoiser),mean(psnr_gaptv),mean(ssim_gaptv),tgaptv);

%% [2.2] DeSCI (with GAP-TV for initialization), TPAMI'18
para.denoiser = 'wnnm'; % WNNM denoising
  para.wnnm_int_fwise = true; % enable GAP-WNNM integrated (with frame-wise denoising)
    para.blockmatch_period = 20; % period of block matching
  para.sigma   = [12]/MAXB; % noise deviation (to be estimated and adapted)
  para.vrange  = 1; % value range
  para.maxiter = [40];
  para.patchsize = 32; % patch size
  para.iternum = 1; % iteration number in WNNM
  para.enparfor = true; % enable parfor
  if para.enparfor % if parfor is enabled, start parpool in advance
      delete(gcp('nocreate')); % delete current parpool
      mycluster = parcluster('local');
      div = 1;
      while nmask/div > mycluster.NumWorkers
          div = div+1;
      end
      parnum = max(min(ceil(nmask/div),mycluster.NumWorkers),1); 
      parpool(mycluster,parnum);
  end

[vdesci,psnr_desci,ssim_desci,tdesci,psnrall] = ...
    gapdenoise_cacti(mask,meas,orig,vgaptv,para); % vgaptv as initialization
  
  tdesci = tdesci + tgaptv;
  delete(gcp('nocreate')); % delete parpool

fprintf('DeSCI mean PSNR %2.2f dB, mean SSIM %.4f, total time % 4.1f s.\n',...
    mean(psnr_desci),mean(ssim_desci),tdesci);

%% [3] save results as mat file
matdir = [resultdir '/savedmat'];
if ~exist(matdir,'dir')
    mkdir(matdir);
end
save([matdir '/desci_' para.type '_' para.name '.mat']);
