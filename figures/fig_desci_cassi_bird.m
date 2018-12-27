%FIG_DESCI_CASSI_BIRD Demonstrate decompress snapshot compressive imaging 
%(DeSCI) results of simulated coded aperture snapshot spectral imaging 
%(CASSI) `bird` dataset.
% Note: Run this demonstration script after completing the reconstruction
%       process TEST_DESCI_CASSI_BIRD and obtaining the saved .mat file.
%   See also TEST_DESCI_CASSI_BIRD, GAPDENOISE_CACTI, GAPDENOISE.
clear; clc;
% close all
addpath('../'); % repository root
addpath(genpath('../packages/')); % packages
% [1] run the corresponding test file for recovery or load the saved
% results
% test_desci_cassi_bird; % run the corresponding test file for recovery
load('../results/savedmat/desci_cassi_bird.mat'); % load the saved results

methods   = {'gaptv','desci'}; % all methods for comparison
methnames = {'GAP-TV','DeSCI'}; % corresponding names of the methods
nim = nframe*nmask; % number of images

%% [2.1] demonstrate the reconstructed spectum
colors = {'Yellow bird', 'Green bird', 'Blue bird', 'Purple bird'}; % colors of the four birds
crops  =  [ 79 413 60 60;  % yellow [x y width height] as insertShape
           356 500 60 60;  % green
           780 460 60 60;  % blue
           950 600 60 60]; % purple
locos = {'southeast','southeast','southeast','northwest'};
       
markers   = {'bx--','r*-'}; % markers

linewidth = 2;
markersize = 4;
textfontsize = 14;
legendfontsize = 12;

% Change default axes fonts.
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultAxesFontSize',textfontsize);
set(0,'DefaultAxesFontWeight','normal');
% Change default text fonts.
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontSize',textfontsize);
set(0,'DefaultTextFontWeight','normal');

% legends
lgall = 'lgtruth';
lgtruth = 'Ground truth';
for imeth = 1:length(methods)
    meth = methods{imeth}; % method
    eval(sprintf('lg%s=''%s'';',meth,methnames{imeth}));
    lgall = [lgall ',lg' meth];
end

f = figure('position',[25 50 800 900]); % PSNR plot
h = tight_subplot(3,2,[0.07 0.08],[.07 .035],[.1 .02]);
axes(h(1));
imshow(insertShape(orig_rgb,'Rectangle',crops,'Color','cyan','Linewidth',8));
title('Original');
axes(h(2));
imshow(meas,[]);
title('Coded frame');
for ic = 1:length(colors)
    color = colors{ic};
    crop = crops(ic,:);
    sp_truth = squeeze(sum(sum(orig(crop(2):crop(2)+crop(4)-1,crop(1):crop(1)+crop(3)-1,:)/MAXB,2),1));
    axes(h(ic+2));
    plot(wavelength,sp_truth/max(sp_truth),'k-','linewidth',linewidth,'markersize',markersize); hold on;
    
    for imeth = 1:length(methods)
        meth = methods{imeth};
        eval(sprintf('sp_%s = squeeze(sum(sum(v%s(crop(2):crop(2)+crop(4)-1,crop(1):crop(1)+crop(3)-1,:),2),1));',meth,meth));
        eval(sprintf('plot(wavelength,rescale(sp_%s,min(sp_truth)/max(sp_truth),1),markers{imeth},''linewidth'',linewidth,''markersize'',markersize);',meth));
        eval(sprintf('corrall(ic,imeth)=corr(sp_%s,sp_truth);',meth));
        % legend with correlation
        eval(sprintf('lg%s=''%s, corr: %.5f'';',meth,methnames{imeth},corrall(ic,imeth)));
    end
    xlim([390 700]); xticks(400:100:700);
    title(color);
    if ic==3 || ic==4
        xlabel('Wavelength (nm)');
    end
    if ic==1 || ic==3
        ylabel('Intensity (a.u.)');
    end
    eval(sprintf('hlg=legend(%s);',lgall));
    set(hlg,'FontSize',legendfontsize);
    legend('location',locos{ic}); 
    legend('boxoff');
end

%% [2.2] show exemplar frames (wavelength)
bands = [9 11 16 21];

nb = length(bands);
nmeth = length(methods);
f = figure('Position',[50 50 1010 950]);
h = tight_subplot(nb,nmeth+1,[.002 .002],[.03 .03],[.02 .02]);
set(f, 'Color', 'white');
for ib = 1:nb
    band = bands(ib);
    lambda = wavelength(band);
    wlmap = (gray*kron(ones(3,1),spectrumRGB(lambda))); % colormap with the corresponding RGB wavelength
    wlmap = wlmap/max(wlmap(:));
    
    axes(h(1+(ib-1)*(nmeth+1)));
    imshow(orig(:,:,band)/MAXB,'colormap',wlmap);
    if ib==1
        title('Ground truth');
    end
    for imeth = 1:nmeth
        meth = methods{imeth};
        axes(h(imeth+1+(ib-1)*(nmeth+1)));
        eval(sprintf('imshow(v%s(:,:,band),''colormap'',wlmap);',meth));
        if ib==1
            title(sprintf('%s',methnames{imeth}));
        end
    end
end
        
%% [2.3] show all the frames wavelength-by-wavelength
%         via vshowSpectralData
addpath(genpath('../packages/'));

opts = [];
  opts.npcol      = 6; % number of columns for subplots
  % opts.magsize    = 0.7; % magnification of the image size
  opts.width      = 1080; % width of the plot
  opts.textpos    = [600 600]; % position of the text on the sub-images
  opts.fontname   = 'Myriad Pro'; % font name of the text on sub-images
  opts.fontsize   = 15; % fone size of the text on sub-images
  opts.fontweight = 'bold'; % font weight of the text on sub-images
% show ground truth
vshowSpectralData(orig/MAXB,wavelength,opts);
for imeth = 1:length(methods)
    meth = methods{imeth};
    eval(sprintf('vshowSpectralData(v%s,wavelength,opts);',meth));
end
