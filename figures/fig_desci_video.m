%DEMO_DESCI_VIDEO Video demonstration of DeSCI reconstruction results.
%   See also TEST_DESCI, GAPDENOISE_CACTI, GAPDENOISE.
clear; clc;
% close all;
% [1] run the test_desci to get the results or load the saved results
% test_desci; % run the test_desci to get the results
load('../results/savedmat/desci_cacti_kobe.mat'); % load the saved results

nim = nframe*nmask; % number of images
methods   = {'gaptv','desci'}; % all methods for comparison
% corresponding names of the methods
methnames = {'GAP-TV','DeSCI'}; % names

% [2] video comparison (.mp4 video and .gif animated image)
framerate = 5;

videodir = '../results/video';
if ~exist(videodir,'dir')
    mkdir(videodir);
end
vname = sprintf('%s/desci_gaptv_%s%d.mp4',videodir,para.name,nim);
gifname = sprintf('%s/desci_gaptv_%s%d.gif',videodir,para.name,nim);
vobj = VideoWriter(vname,'MPEG-4');
vobj.FrameRate = framerate;
open(vobj);
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);
f = figure('Position',[50 100 800 800]);
set(f, 'Color', 'white');
for kv = 1:nim
    subplot(221); % (a) coded frame
    imshow(uint8(imnorm(meas(:,:,floor((kv-1)/nmask)+1))*MAXB));
    title(sprintf('\\rm Coded frame #\\bf%d',floor((kv-1)/nmask)+1));
    subplot(222); % (b) original video
    imshow(uint8(orig(:,:,kv)));
    title(sprintf('\\rm Original #%d',kv));
    for imeth = 1:length(methods)
        % subplot
        meth = methods{imeth};
        subplot(2,2,imeth+2);
        eval(sprintf('imshow(uint8(v%s(:,:,kv)*MAXB));',meth));
        eval(sprintf('impsnr=psnr_%s(kv);imssim=ssim_%s(kv);',meth,meth));
        if imeth == length(methods)
            title(sprintf('\\rm%s (\\bf%2.2f\\rm dB, \\bf%.4f\\rm)',...
                methnames{imeth},impsnr,imssim));
        else
            title(sprintf('\\rm%s (%2.2f dB, %.4f)',...
                methnames{imeth},impsnr,imssim));
        end
    end
    
    frame = getframe(f);
    writeVideo(vobj,frame);
    [vgif,map] = rgb2ind(frame.cdata,256,'nodither');
    if kv==1 % save as gif (first frame)
        imwrite(vgif,map,gifname,'DelayTime',1/framerate,'LoopCount',inf); 
    else % save as gif
        imwrite(vgif,map,gifname,'DelayTime',1/framerate,...
            'WriteMode','append');
    end
end
close(vobj);

