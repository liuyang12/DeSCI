function [  ] = vshowSpectralData( vdata, wavelength, opts )
%VSHOWSPECTRALDATA Showing a set of spectral data volume as suplots of a
%figure.
%   See also SPECTRUMRGB, TIGHT_SUBPLOT.

% default parameter setting
npcol    = 6; % number of columns for subplots
% magsize  = 1; % magnification of the image size
width    = 1080; % width of the plot
textpos  = [12 24]; % position of the text on the sub-images
fontsize = 15; % font size of the text (wavelength nm) on images
fontweight = 'normal'; % font weight of the text (wavelength) on images
fontname = 'FixedWidth'; % font name of the text (wavelength nm) on images
% parameters from the input opts
if isfield(opts,'npcol'),           npcol = opts.npcol;      end
if isfield(opts,'magsize'),       magsize = opts.magsize;    end
if isfield(opts,'width'),           width = opts.width;      end
if isfield(opts,'textpos'),       textpos = opts.textpos;    end
if isfield(opts,'fontsize'),     fontsize = opts.fontsize;   end
if isfield(opts,'fontweight'), fontweight = opts.fontweight; end
if isfield(opts,'fontname'),     fontname = opts.fontname;   end

% nw = length(wavelength); % number of wavelength bands
[nr,nc,nw] = size(vdata); % number of rows, columns, and wavelength bands
nprow = ceil(nw/npcol); % number of rows for subplots
% f = figure;
% f = figure('Position',[50 50 1100 900]); % depending on the size of images
% f = figure('Position',[50 50 npcol*nc*magsize nprow*nr*magsize]); % depending on the size of images
f = figure('Position',[50 50 width width/(npcol*nc)*nprow*nr]); % depending on the size of images
h = tight_subplot(nprow,npcol,[.002 .002],[.01 .01],[.01 .01]);
set(f,'Color','white');
% display each wavelength in each black with the same order as the
% wavelength
for iw = 1:nw
    lambda = wavelength(iw); % wavelength (nm)
    % colormap with the corresponding RGB wavelength
    wmap = (gray*kron(ones(3,1),spectrumRGB(lambda)));
    wmap = wmap/max(wmap(:));
    
    % subplot(nprow,npcol,iw);
    axes(h(iw));
    imshow(vdata(:,:,iw),'Colormap',wmap);
    text(textpos(1),textpos(2),[num2str(lambda,3) ' nm'],'Color','white',...
        'FontWeight',fontweight,'FontName',fontname,'FontSize',fontsize);
end
% display white image for the remaining blacks
for iw = nw+1:nprow*npcol
    axes(h(iw));
    imshow(ones(size(vdata,2)));
end

end

