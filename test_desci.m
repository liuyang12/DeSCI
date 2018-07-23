%TEST_DESCI Test decompress snapshot compressive imaging (DeSCI) for 
%simulated coded aperture compressive temporal imaging (CACTI) dataset.
% Reference
%   [1] Y. Liu, X. Yuan, J. Suo, D.J. Brady, and Q. Dai, Rank Minimization
%       for Snapshot Compressive Imaging, preprint, 2018.
%   [2] X. Yuan, Generalized alternating projection based total variation 
%       minimization for compressive sensing, in ICIP, pp. 2539-2543, 2016.
% Contact
%     Xin Yuan, Bell Labs, xyuan@bell-labs.com, initial version Jul 2, 
%       2015.
%     Yang Liu, Tsinghua University, y-liu16@mails.tsinghua.edu.cn, last
%       update Jul 13, 2018.
%   See also GAPDENOISE.
clear; clc;
% close all
% [0] environment configuration
addpath(genpath('./algorithms')); % algorithms
addpath(genpath('./packages')); % packages
addpath(genpath('./utils')); % utilities

datasetdir = './dataset'; % dataset
resultdir  = './results'; % results

% [1] load dataset
para.type   = 'cacti'; % type of dataset, cassi or cacti
para.name   = 'kobe'; % name of dataset
para.number = 32; % number of frames in the dataset

datapath = sprintf('%s/%s%d_%s.mat',datasetdir,para.name,...
    para.number,para.type);

load(datapath); % mask, meas, orig (and para)

para.nframe = 1; % number of coded frames in this test
para.MAXB   = 255;

[nrow,ncol,nmask] = size(mask);
nframe = para.nframe; % number of coded frames in this test
MAXB = para.MAXB;

% [1.2] parameter setting for GAP-TV and GAP-WNNM
para.Mfunc  = @(z) A_xy(z,mask);
para.Mtfunc = @(z) At_xy_nonorm(z,mask);

para.Phisum = sum(mask.^2,3);
para.Phisum(para.Phisum==0) = 1;

para.enparfor = true; % enable parfor for multi-CPU acceleration
para.flag_iqa = false; % disable image quality assessments in iterations

%% [2.1] GAP-TV(-acc), ICIP2016
para.lambda   =    1; % correction coefficiency
para.maxiter  =  100; % maximum iteration
para.acc      =    1; % enable acceleration
para.denoiser = 'tv'; % TV denoising
  para.tvweight = 0.07*255/MAXB; % weight for TV denoising
  para.tviter   = 5; % number of iteration for TV denoising
  
[vgaptv,psnr_gaptv,ssim_gaptv,tgaptv] = ...
    gapdenoise_cacti(mask,meas,orig,[],para);

fprintf('GAP-%s mean PSNR %2.2f dB, mean SSIM %.4f, total time % 4.1f s.\n',...
    upper(para.denoiser),mean(psnr_gaptv),mean(ssim_gaptv),tgaptv);

%% [2.2] DeSCI, 2018
para.acc = 1; % enable acceleration
para.denoiser = 'wnnm'; % WNNM denoising
  para.wnnm_int = true; % enable GAP-WNNM integrated
    para.blockmatch_period = 20; % period of block matching
  para.sigma   = [100 50 25 12   6]/MAXB; % noise deviation (to be estimated and adapted)
  para.vrange  = 1; % value range
  para.maxiter = [ 60 60 60 60  60];
  para.iternum = 1; % iteration number in WNNM
  para.enparfor = true; % enable parfor
  if para.enparfor % if parfor is enabled, start parpool in advance
      mycluster = parcluster('local');
      delete(gcp('nocreate')); % delete current parpool
      poolobj = parpool(mycluster,min(nmask,mycluster.NumWorkers));
  end

[vdesci,psnr_desci,ssim_desci,tdesci,psnrall] = ...
    gapdenoise_cacti(mask,meas,orig,[],para);
  
  delete(poolobj); % delete pool object

fprintf('DeSCI mean PSNR %2.2f dB, mean SSIM %.4f, total time % 4.1f s.\n',...
    mean(psnr_desci),mean(ssim_desci),tdesci);

%% [3] save results as mat file
matdir = [resultdir '/savedmat'];
if ~exist(matdir,'dir')
    mkdir(matdir);
end
save([matdir '/desci_' para.name num2str(nframe*nmask) '.mat']);

  