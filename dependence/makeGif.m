function makeGif(volume, delayTime, filename)
%%% Function to make a gif from a 3D volume which can also have colour info (4 dims)
%%% or no colour info (3dims) - P. Pandey
close hidden all
if ndims(volume) == 4
    n = size(volume,4);
elseif ndims(volume) == 3
    n = size(volume,3);
else
    error('Incorrectly sized input volume')
    return
end
for idx = 1:n
    if ndims(volume) == 4
        imshow(volume(:,:,:,idx));
    else
        imshow(volume(:,:,idx));
    end
    drawnow
    frame = getframe(1);
    im{idx} = frame2im(frame);
end
for idx = 1:n
    [gif,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(gif,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(gif,map,filename,'gif','WriteMode','append','DelayTime',delayTime);
    end
end
close all