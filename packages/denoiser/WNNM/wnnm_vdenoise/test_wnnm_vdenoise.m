%TEST_WNNM_VDENOISE Test the weighted nuclear norm minimization 
%(WNNM)-based video denoising algorithm.
%   See also WNNM_VDENOISE, VPARACONFIG.
clear; clc;
% close all
% [0] environment configuration
addpath(genpath('./packages')); % packages
addpath('./utils'); % utilities

vdir = './video'; % video directory
rdir = './result'; % result directory
  rvdir  = [rdir '/video']; % video result directory
  rimdir = [rdir '/image']; % image result directory

% [1] read original video and add Gaussian noise to the video
vname      = 'bicycle'; % name of the video
GRAYSCALE  = true; % grayscale video [forced]
sigma      = 20; % noise deviation w.r.t. 8-bit grayscale (0~255)
framerate  = 5; % framerate for saving videos
fnum       = [1:8]; % number of the frames used in the test
cfno       = round(length(fnum)/2); % center frame number

videopath = sprintf('%s/%s.avi',vdir,vname);

orgpath = sprintf('%s/%s_org%d-%d.avi',rvdir,vname,fnum(1),fnum(end));
nospath = sprintf('%s/%s_nos%d-%d_sigma%d.avi',rvdir,vname,fnum(1),...
    fnum(end),sigma);

if exist(orgpath,'file') && exist(nospath,'file') % load exist file
    orgv = readVideo(orgpath);
    nosv = readVideo(nospath);
else % re-generate noisy data
    originalvideo = uint8(readVideo(videopath));
    otherdims = repmat({':'},1,ndims(originalvideo)-1);
    if GRAYSCALE && ndims(originalvideo)>=4
        orgv = vrgb2gray(originalvideo(otherdims{:},fnum));
    else
        orgv = originalvideo(otherdims{:},fnum);
    end

    nosv = sigma*randn(size(orgv))+double(orgv); % add Gaussian noise
    % save the denoised video to .avi file
    if ~exist(rdir,'dir')
        mkdir(rdir);
    end
    if ~exist(rvdir,'dir')
        mkdir(rvdir);
    end
    write_video(orgv,framerate,orgpath);
    write_video(nosv,framerate,nospath);
end

% [2] apply WNNM-based video denoising
orgv = double(orgv);
%% [2.0] noisy indexes
[meanpsnr_noisy,psnr_noisy] = vpsnr(nosv,orgv,255);
[meanssim_noisy,ssim_noisy] = vssim(nosv,orgv);

%% [2.1] V-BM3D denoising
tic
[~,vbm3dv] = VBM3D(nosv/255,sigma,0,0); % with input noise deviation
tvbm3d = toc;
vbm3dv = vbm3dv*255;
[meanpsnr_vbm3d,psnr_vbm3d] = vpsnr(double(vbm3dv),orgv,255);
[meanssim_vbm3d,ssim_vbm3d] = vssim(double(vbm3dv),orgv);

%% [2.2] V-BM4D denoising
tic
vbm4dv = vbm4d(nosv,sigma,'lc',1,1,1,0); % with input standard deviation
% vbm4dv = vbm4d(nosv,-1,'lc',1,1,1,0); % enable noisy estimation
tvbm4d = toc;
[meanpsnr_vbm4d,psnr_vbm4d] = vpsnr(double(vbm4dv),orgv,255);
[meanssim_vbm4d,ssim_vbm4d] = vssim(double(vbm4dv),orgv);

%% [2.3] WNNM-based video denoising (C-style implementation, structure)
vpara.sigma = sigma*2; % sigma*2 optimal for bus-8
vpara.enparfor = true; % enable parfor
tic
wnnm_cv = wnnmvdenoiser(nosv,[],[],vpara);
twnnm_c = toc;
[meanpsnr_wnnm_c,psnr_wnnm_c] = vpsnr(wnnm_cv,orgv,255);
[meanssim_wnnm_c,ssim_wnnm_c] = vssim(wnnm_cv,orgv);

%% [2.4] WNNM-based video denoising (MATLAB-style implementation, matrix)
vpara = vparaconfig(sigma);
  vpara.adaptboost = false; % enable adaptice boosting (AB) in WNNM
  vpara.abeta      = 1; % eta in AB noise estimation
  
tic
wnnmv = wnnm_vdenoise_par(nosv,orgv,vpara);
twnnm = toc;
[meanpsnr_wnnm,psnr_wnnm] = vpsnr(wnnmv,orgv,255);
[meanssim_wnnm,ssim_wnnm] = vssim(wnnmv,orgv);

%% [3] save the result as .mat file
mdir = [rdir '/savedmat']; % saved .mat file directory
save([mdir '/wnnmv_' vname num2str(length(fnum)) '_'...
    char(datetime('today','format','yyMMdd')) '.mat']);

%% [4] result demonstration
% [4.1] figure
linewidth = 1.5;
markersize = 6;
figure; % PSNR plot
plot(psnr_vbm3d,'ko--','linewidth',linewidth,'markersize',markersize); 
hold on; % grid on;
plot(psnr_vbm4d,'k<--','linewidth',linewidth,'markersize',markersize);
plot(psnr_wnnm_c,'kd-','linewidth',linewidth,'markersize',markersize);
plot(psnr_wnnm,'ks-','linewidth',linewidth,'markersize',markersize);
title(sprintf('PSNR (dB) of %s-%d',vname,length(fnum)));
set(gcf,'position',[600 500 440 400]);
xlim([1 length(fnum)]);
figure; % SSIM plot
plot(ssim_vbm3d,'ko--','linewidth',linewidth,'markersize',markersize); 
hold on; % grid on;
plot(ssim_vbm4d,'k<--','linewidth',linewidth,'markersize',markersize);
plot(ssim_wnnm_c,'kd-','linewidth',linewidth,'markersize',markersize);
plot(ssim_wnnm,'ks-','linewidth',linewidth,'markersize',markersize);
title(sprintf('SSIM of %s-%d',vname,length(fnum)));
set(gcf,'position',[600 500 700 400]);
xlim([1 length(fnum)]);
lgvbm3d = sprintf('V-BM3D (%2.2f dB, %.4f)',mean(psnr_vbm3d),mean(ssim_vbm3d));
lgvbm4d = sprintf('V-BM4D (%2.2f dB, %.4f)',mean(psnr_vbm4d),mean(ssim_vbm4d));
lgwnnm_c = sprintf('WNNM-C (%2.2f dB, %.4f)',mean(psnr_wnnm_c),mean(ssim_wnnm_c));
lgwnnm = sprintf('WNNM (%2.2f dB, %.4f)',mean(psnr_wnnm),mean(ssim_wnnm));
hlg = legend(lgvbm3d,lgvbm4d,lgwnnm_c,lgwnnm); legend('boxoff');
set(hlg,'FontSize',11);
legend('location','northeastoutside');
legend('boxoff');

%% [4.2] save all videos in a single window both in mp4 and gif
if ~exist(rdir,'dir')
    mkdir(rdir);
end
if ~exist(rvdir,'dir')
    mkdir(rvdir);
end
framerate = 5;
mp4name = sprintf('%s/%s%d_wnnm_%s.mp4',rvdir,vname,length(fnum),...
    char(datetime('today','format','yyMMdd')));
gifname = sprintf('%s/%s%d_wnnm_%s.gif',rvdir,vname,length(fnum),...
    char(datetime('today','format','yyMMdd')));
vobj = VideoWriter(mp4name,'MPEG-4');
vobj.FrameRate = framerate;
open(vobj); % assert to close the vobj after writing the video
f = figure('Position',[50 100 900 680]);
h = tight_subplot(2,3,[.01 .01],[.02 .02],[.02 .02]);
set(f, 'Color', 'white');
for k = 1:length(fnum)
    axes(h(1));
    imshow(uint8(orgv(:,:,k)));
    title(sprintf('(a) Original #%d',fnum(k)));
    axes(h(2));
    imshow(uint8(nosv(:,:,k)));
    title(sprintf('(b) Noisy \\sigma = %d',sigma));
    axes(h(3));
    imshow(uint8(vbm3dv(:,:,k)));
    title(sprintf('(c) V-BM3D, PSNR %2.2f dB',psnr_vbm3d(k)));
    axes(h(4));
    imshow(uint8(vbm4dv(:,:,k)));
    title(sprintf('(d) V-BM4D, PSNR %2.2f dB',psnr_vbm4d(k)));
    axes(h(5));
    imshow(uint8(wnnm_cv(:,:,k)));
    title(sprintf('(e) WNNM-C, PSNR %2.2f dB',psnr_wnnm_c(k)));
    axes(h(6));
    imshow(uint8(wnnmv(:,:,k)));
    title(sprintf('(f) WNNM, PSNR %2.2f dB',psnr_wnnm(k)));
    
    frame = getframe(f);
    writeVideo(vobj,frame);
    [vgif,map] = rgb2ind(frame.cdata,256,'nodither');
    if k==1
        imwrite(vgif,map,gifname,'DelayTime',1/framerate,'LoopCount',inf); % save as gif
    else
        imwrite(vgif,map,gifname,'DelayTime',1/framerate,'WriteMode','append'); % save as gif
    end
end
close(vobj);


