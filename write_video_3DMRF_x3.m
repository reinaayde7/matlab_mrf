function write_video_3DMRF_x3(Image,filename,framerate,colormapP, climst)
% write_gif creates a gif.
% saves it in filename.
% time between frames = timer.
% BAM !
% Image is a [ymax,xmax,nt,sets] matrix.
% filename ex: 'My_animation.gif'.
% papersize = [vertical horizontal];

if isempty(colormapP)
    colormapP = 'gray';
end

fig=figure; set(gcf,'color','w');
fig.Position = [0 0 2000 800];
set(gcf,'windowstyle','normal');


[~,~,nt,sets] = size(Image);

% % % % % % --------------------- VIDEO GENERATOR -----------------------------
% % % % % writerObj=VideoWriter(filename); % Change the video file name
% % % % % %Change the frame rate as per requirements
% % % % % writerObj.FrameRate = framerate; 
% % % % % open(writerObj);
% % % % % 
% % % % %     for n = 1:nt
% % % % %         t = tiledlayout(1, 3, "TileSpacing", "none");
% % % % %         for s=1:sets
% % % % %             nexttile;
% % % % %             ImageAux = Image(:,:,n,s);
% % % % %             imagesc(ImageAux, climst), colormap(colormapP)
% % % % %             axis off
% % % % %         end
% % % % % 
% % % % %         h = axes(fig,'visible','off'); 
% % % % %         c = colorbar(h,'Position',[0.94 0.168 0.022 0.7], 'FontSize',20);  % attach colorbar to h
% % % % %         colormap(c,colormapP)
% % % % %         caxis(h,climst);  
% % % % % 
% % % % %         frame = getframe(fig);
% % % % %         writeVideo(writerObj,frame);
% % % % % 
% % % % %     end
% % % % % 
% % % % % writerObj=VideoWriter(filename); % Change the video file name
% % % % % %Change the frame rate as per requirements
% % % % % writerObj.FrameRate = framerate; 
% % % % % open(writerObj);

% --------------------- GIF GENERATOR -----------------------------
    for n = 1:nt
        t = tiledlayout(1, 3, "TileSpacing", "none");
        for s=1:sets
            ax = subplot(1,sets,s);
            nexttile;
            ImageAux = Image(:,:,n,s);
            imagesc(ImageAux, climst), colormap(ax, colormapP)
            imagesc(ImageAux, climst), colormap(colormapP)
            axis off
        end
        
        h = axes(fig,'visible','off'); 
        c = colorbar(h,'Position',[0.94 0.168 0.022 0.7], 'FontSize',20);  % attach colorbar to h
        colormap(c,colormapP)
        caxis(h,climst);  
    
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',timer);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',timer);
        end

        writeVideo(writerObj,frame);

    end




end
