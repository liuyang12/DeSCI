%FIG_DESCI_CASSI Demonstrate decompress snapshot compressive imaging 
%(DeSCI) results of simulated coded aperture snapshot spectral imaging 
%(CASSI) `toy` dataset.
% Note: Run this demonstration script after completing the reconstruction
%       process TEST_DESCI_CASSI and obtaining the saved .mat file.
%   See also TEST_DESCI_CASSI, GAPDENOISE_CACTI, GAPDENOISE.
clear; clc;
% close all
addpath('../'); % repository root
addpath(genpath('../packages/')); % packages
% [1] run the corresponding test file for recovery or load the saved
% results
% test_desci_cassi % run the corresponding test file for recovery
load('../results/savedmat/desci_cassi_toy.mat'); % load the saved results

methods   = {'gaptv','desci'}; % all methods for comparison
methnames = {'GAP-TV','DeSCI'}; % corresponding names of the methods
nim = nframe*nmask; % number of images

%% [2.1] demonstrate the reconstructed spectum
colors = {'Brown', 'Orange', 'Blue', 'Red'}; % colors of the four birds
crops  =  [150 280 24 24;  % brown [x y width height] as insertShape
           340 210 24 24;  % orange
           349 396 24 24;  % blue
           349 439 24 24]; % red
locos = {'northwest','northwest','northwest','northwest'};
       
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

f = figure('position',[25 50 700 900]); % PSNR plot
h = tight_subplot(3,2,[0.07 0.08],[.07 .035],[.1 .02]);
% subplot(3,2,1); % show crops on original RGB image
axes(h(1));
imshow(insertShape(orig_rgb,'Rectangle',crops,'Color','cyan','Linewidth',5));
title('Original');
% subplot(3,2,2); % coded measurement
axes(h(2));
imshow(meas,[]);
title('Coded frame');
for ic = 1:length(colors)
    color = colors{ic};
    crop = crops(ic,:);
    sp_truth = squeeze(sum(sum(orig(crop(2):crop(2)+crop(4)-1,crop(1):crop(1)+crop(3)-1,:)/MAXB,2),1));
    % figure;
    % subplot(3,2,ic+2);
    axes(h(ic+2));
    plot(wavelength,sp_truth/max(sp_truth),'k-','linewidth',linewidth,'markersize',markersize,'DisplayName','Ground truth'); hold on;
    
    for imeth = 1:length(methods)
        meth = methods{imeth};
        eval(sprintf('sp_%s = squeeze(sum(sum(v%s(crop(2):crop(2)+crop(4)-1,crop(1):crop(1)+crop(3)-1,:),2),1));',meth,meth));
        eval(sprintf('corrall(ic,imeth)=corr(sp_%s,sp_truth);',meth)); % corr is scale invariant 
        eval(sprintf('plot(wavelength,rescale(sp_%s,min(sp_truth)/max(sp_truth),1),markers{imeth},''linewidth'',linewidth,''markersize'',markersize,''DisplayName'',[methnames{imeth} '', corr: %s'']);',meth,num2str(corrall(ic,imeth),'%.4f')));
    end
    title(color);
    if ic==3 || ic==4
        xlabel('Wavelength (nm)');
    end
    if ic==1 || ic==3
        ylabel('Intensity (a.u.)');
    end
    ylim([0 1]);
    hlg = legend('show');
    set(hlg,'FontSize',legendfontsize);
    legend('location',locos{ic}); 
    legend('boxoff');
end

%% [2.2] show exemplar frames (wavelength)
bands = [1 10 20 31];

nb = length(bands);
nmeth = length(methods);
f = figure('Position',[50 50 700 950]);
h = tight_subplot(nb,nmeth+1,[.002 .002],[.03 .03],[.02 .02]);
set(f, 'Color', 'white');
for ib = 1:nb
    band = bands(ib);
    lambda = wavelength(band);
    wlmap = (gray*kron(ones(3,1),spectrumRGB(lambda))); % colormap with the corresponding RGB wavelength
    wlmap = wlmap/max(wlmap(:));
    
    axes(h(1+(ib-1)*(nmeth+1)));
    imshow(uint8(orig(:,:,band)),'colormap',wlmap);
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

