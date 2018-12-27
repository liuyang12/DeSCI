%TEST_WNNM_IMDENOISE Test the weighted nuclear norm minimization 
%(WNNM)-based image denoising function.
%   See also WNNM_IMDENOISE, PARACONFIG.
clear; clc;
% close all
addpath(genpath('./apg_partial'));
% [1] parameter configuration
nsigma = 50; % standard deviation of the noise

% [2] read original image (and noisy image, possibly)
orgim = double(imread('boat.png'));
% % nosim = double(imread('boat_n100.png'));
% orgim = double(rgb2gray(imread('bus25.png')));
% nosim = double(rgb2gray(imread('bus25_pref.png')));
nosim = orgim + nsigma*randn(size(orgim)); % add white Gaussian noise

psnr0 = psnr(nosim,orgim,255);
fprintf('Noisy image with nsigma = %2.1f, PSNR = %2.2f dB.\n',nsigma,psnr0);

% para.iternum = 1;
% [3] perform WNNM-based image denoising
para = paraconfig(nsigma); % parameter configuration according to noise level
  para.adaptboost = false;  % enable adaptive boosting (AB) for WNNM
  para.abeta      = 1;     % eta for noise estimation in WNNM-AB
tic
estim = wnnm_imdenoise(nosim,orgim,para); toc

psnr1 = psnr(estim,orgim,255);
fprintf('WNNM denoised image with nsigma = %2.1f, PSNR = %2.2f dB.\n',nsigma,psnr1);

%% [4] visulization
figure;
subplot(221); imshow(orgim,[]); title('Original image (bus25)');
subplot(222); imshow(nosim,[]); title(sprintf('Noisy image (\\sigma=%2.1f, PSNR=%2.2f dB)',nsigma,psnr0));
subplot(224); imshow(estim,[]); title(sprintf('WNNM image (PSNR= %2.2f dB)',psnr1));
