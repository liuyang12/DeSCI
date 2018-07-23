function y = read_video(file_name, frames)
    fprintf('Loading test sequence "%s" \n', file_name);
    mov       = VideoReader(file_name);
    if ~exist('frames','var') || isempty(frames)
        frames = mov.NumberOfFrames;
    end
    nFrames   = min(frames,mov.NumberOfFrames);
    vidHeight = mov.Height;
    vidWidth  = mov.Width;
    isRGB = strcmpi(mov.VideoFormat,'RGB24');
    if isRGB
        y = zeros(vidHeight,vidWidth,3,nFrames);
    else
        y = zeros(vidHeight,vidWidth,nFrames);
    end
    for n = 1:nFrames,
        frame = read(mov, n);
        if isRGB
            y(:,:,:,n)  = frame;
        else
            y(:,:,n)  = frame;
        end
    end
end