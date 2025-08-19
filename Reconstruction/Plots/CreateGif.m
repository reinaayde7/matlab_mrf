% Script to generate gifs from 3D matrices
% Jesus Fajardo jesuserf@umich.edu

filename = strcat('PM_SVDcompressed_Nex1000','.mat'); %mat filename
load(filename,'t1map3d','t2map3d','m0map3d');
t1map3d = rot90((reshape(t1map3d,[size(t1map3d,1) size(t1map3d,1) 1 size(t1map3d,3)])));
t2map3d = rot90((reshape(t2map3d,[size(t2map3d,1) size(t2map3d,1) 1 size(t2map3d,3)])));
m0map3d = rot90((abs(reshape(m0map3d,[size(m0map3d,1) size(m0map3d,1) 1 size(m0map3d,3)]))));

addpath('./Colormap');
addpath('./dependence');

D = (t2pf(:,:,:));
%montage(D, "Indices",40:45,"DisplayRange",[0 300], 'hot');
%montage(D, "Indices",29:34,"DisplayRange",  [0, 300]);
climst1 = [0, 3000];
climst2 = [0, 300];
% Set up the 3D array and define the axis and speed
axis = 3;
speed = 5;

% Load the 3D matrix from .mat file
figure(1)
filename = 'AxT1.gif';
check_flag = 0;
for n = 10:50
    subplot();%'Position', [0., 0., 1., 1.]
    data = squeeze(D(160:310,140:290,n));
    imagesc(data,climst2);%axis image;colormap gray;
    axis off;
    colormap(T2colormap);
    %title(strcat(num2str(ii),'/',num2str(NumFrames), ' slice'));
    daspect([1 1 1])
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      
    if check_flag == 0
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
        check_flag = 1;
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end

end
