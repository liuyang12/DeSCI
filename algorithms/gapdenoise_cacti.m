function [v_,psnr_,ssim_,t_,psnrall] = gapdenoise_cacti(mask,meas,orig,v0,para)
%GAPDENOISE_CACTI GAP-Denoise frame for recontruction of  CACTI high-speed 
%imaging.
%   See also GAPDENOISE, TEST_GAPDENOISE.
iframe   = 1; % start frame number
maskdirection  = 'plain'; % direction of the mask
if isfield(para,'iframe');               iframe = para.iframe; end
if isfield(para,'maskdirection'); maskdirection = para.maskdirection; end
[nrow,ncol,nmask]  = size(mask);
nframe = para.nframe;
MAXB   = para.MAXB;

psnrall = [];
% ssimall = [];
v_ = zeros([nrow ncol nmask*nframe],'like',meas);
tic
% coded-frame-wise denoising
for kf = 1:nframe
    fprintf('GAP-%s Reconstruction frame-block %d of %d ...\n',...
        upper(para.denoiser),kf,nframe);
    if ~isempty(orig)
        para.orig = orig(:,:,(kf-1+iframe-1)*nmask+(1:nmask))/MAXB;
    end
    y = meas(:,:,kf+iframe-1)/MAXB;
    if isempty(v0) % raw initialization
        para.v0 = [];
    else % given initialization
        switch lower(maskdirection)
            case 'plain'
                para.v0 = v0(:,:,(kf-1)*nmask+(1:nmask));
            case 'updown'
                if mod(kf+iframe-1,2) == 0 % even frame (falling of triangular wave)
                    para.v0 = v0(:,:,(kf-1)*nmask+(1:nmask));
                else % odd frame (rising of triangular wave)
                    para.v0 = v0(:,:,(kf-1)*nmask+(nmask:-1:1));
                end
            case 'downup'
                if mod(kf+iframe-1,2) == 1 % odd frame (rising of triangular wave)
                    para.v0 = v0(:,:,(kf-1)*nmask+(1:nmask));
                else % even frame (falling of triangular wave)
                    para.v0 = v0(:,:,(kf-1)*nmask+(nmask:-1:1));
                end
            otherwise
                error('Unsupported mask direction %s!',lower(maskdirection));
        end
        
    end
    if isfield(para,'wnnm_int') && para.wnnm_int % GAP-WNNM integrated
        if isfield(para,'flag_iqa') && ~para.flag_iqa % ImQualAss disabled
            v = gapwnnm_int(y,para);
        else
            [v,psnrall(kf,:)] = gapwnnm_int(y,para);
        end
    elseif isfield(para,'wnnm_int_fwise') && para.wnnm_int_fwise % GAP-WNNM integrated (with frame-wise denoising)
        if isfield(para,'flag_iqa') && ~para.flag_iqa % ImQualAss disabled
            v = gapwnnm_int_fwise(y,para);
        else
            [v,psnrall(kf,:)] = gapwnnm_int_fwise(y,para);
        end
    else
        if isfield(para,'flag_iqa') && ~para.flag_iqa % ImQualAss disabled
            v = gapdenoise(y,para);
        else
            [v,psnrall(kf,:)] = gapdenoise(y,para);
        end
    end
    switch maskdirection
        case 'plain'
            v_(:,:,(kf-1)*nmask+(1:nmask)) = v;
        case 'updown'
            if mod(kf+iframe-1,2) == 0 % even frame (falling of triangular wave)
                v_(:,:,(kf-1)*nmask+(1:nmask)) = v;
            else % odd frame (rising of triangular wave)
                v_(:,:,(kf-1)*nmask+(nmask:-1:1)) = v;
            end
        case 'downup'
            if mod(kf+iframe-1,2) == 1 % odd frame (rising of triangular wave)
                v_(:,:,(kf-1)*nmask+(1:nmask)) = v;
            else % even frame (falling of triangular wave)
                v_(:,:,(kf-1)*nmask+(nmask:-1:1)) = v;
            end
        otherwise 
            error('Unsupported mask direction %s!',lower(maskdirection));
    end
end
t_ = toc;
% image quality assessments
psnr_ = zeros([1 nmask*nframe]);
ssim_ = zeros([1 nmask*nframe]);
if ~isempty(orig)
    for kv = 1:nmask*nframe
        psnr_(kv) = psnr(double(v_(:,:,kv)),double(orig(:,:,kv)/MAXB),max(max(max(double(orig(:,:,kv))/MAXB))));
        ssim_(kv) = ssim(double(v_(:,:,kv)),double(orig(:,:,kv)/MAXB));
    end
end

end

