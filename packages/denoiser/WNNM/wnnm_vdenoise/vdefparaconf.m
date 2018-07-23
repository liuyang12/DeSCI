function [para] = vdefparaconf(para)
%VDEFPARACONF Default parameter configuration for WNNM-based video
%denoising.
%   para = VDEFPARACONF(para) returns the parameters para by adding the
%   fields defined in the default parameters but undefined in the original
%   para.
%   See also VPARACONFIG.
vrange = 255; % default range of the input signal for denoising
if isfield(para,'vrange'), vrange = para.vrange; end
vp = vparaconfig(para.sigma*255/vrange); % default parameters
if isfield(para,'vrange')
    vp.nsigma = para.sigma/vrange;
else
    vp.nsigma = para.sigma;
end
if isfield(para,'patchsize') && ~isfield(para,'patchstep')
    para.patchstep = floor(para.patchsize/2-1);
end
fnames = fieldnames(vp);
for in = 1:length(fnames)
    fn = fnames{in};
    if ~isfield(para,fn)
        eval(sprintf('para.%s=vp.%s;',fn,fn));
    end
end

end

