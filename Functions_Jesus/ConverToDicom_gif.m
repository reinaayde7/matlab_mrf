


filename2 = strcat('PM_SVDcompressed_Nex1000','.mat'); %mat filename
load(filename2,'t1map3d','t2map3d','m0map3d');
t1pf = t1map3d;
t2pf = t2map3d;
t1pf = rot90((reshape(t1pf,[size(t1pf,1) size(t1pf,1) 1 size(t1pf,3)])));
t2pf = rot90((reshape(t2pf,[size(t2pf,1) size(t2pf,1) 1 size(t2pf,3)])));
t1pf = flip(t1pf,4);
t2pf = flip(t2pf,4);

filename = strcat('PM_SVDcompressed_Nex1000','.mat'); %mat filename
load(filename,'t1map3d','t2map3d','m0map3d');
t1map3d = rot90((reshape(t1map3d,[size(t1map3d,1) size(t1map3d,1) 1 size(t1map3d,3)])));
t2map3d = rot90((reshape(t2map3d,[size(t2map3d,1) size(t2map3d,1) 1 size(t2map3d,3)])));
m0map3d = rot90((abs(reshape(m0map3d,[size(m0map3d,1) size(m0map3d,1) 1 size(m0map3d,3)]))));

%{
%% Load Dicom to compare
filename = 'AxT2_DL.dcm';
[info_sample] = dicominfo(filename);
[sample] = dicomread(filename);
sample = squeeze(sample(:,:,1,:));
%}



%file_name1 = strcat('DICOM_T1map','.dcm');
%file_name2 = strcat('DICOM_T2map','.dcm');
%file_name3 = strcat('DICOM_M0map','.dcm');

%dicomwrite(int16(t1map3d), file_name1);%
%dicomwrite(int16(t2map3d), file_name2); %(int16(t2map3d), file_name2, info, 'CreateMode', 'Copy', 'MultiframeSingleFile', 'true')
%dicomwrite(int16(m0map3d), file_name3);

cd('..');
cd('..');
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T1cm.mat');
load('./Colormap/T2cm.mat');
%cd('data/1106')

%m0map3d = flipud(m0map3d);

%hand_drawn_mask_phantom(t1map3d(:,:,7), 14);


%Crop image handrawn
%quant_map = squeeze(m0map3d(90:290,100:300,22));
%{
size_image_x = size(quant_map,1);
size_image_y = size(quant_map,2);
figure, imagesc(1:size_image_x, 1:size_image_y, quant_map);
h = drawfreehand;
ROI = createMask(h);
MaskedImage = ROI.*quant_map;
data = nonzeros(MaskedImage);
%histogram(data,15)
%imshow(MaskedImage,[]);
%hand_drawn_mask_phantom(t2map3d, 1);
%} 

D = (t1pf(:,:,:));
%montage(D, "Indices",40:45,"DisplayRange",[0 300], 'hot');
%montage(D, "Indices",29:34,"DisplayRange",  [0, 300]);
climst1 = [0, 3000];
climst2 = [0, 300];
%{
figure('DefaultAxesFontSize',2)
subplot()
climst1 = [0, 3000];
climst2 = [0, 300];
imagesc(squeeze(D(:,:,38)),climst2);
colormap("gray");
title('T2')
daspect([1 1 1])
colorbar


for i = 230,
figure('DefaultAxesFontSize',2)
subplot()
climst1 = [0, 3000];
climst2 = [0, 300];
imagesc(squeeze(D(130:350,i,10:50)),climst1);
colormap(T1colormap);
title('T1')
daspect([1 1 1])
colorbar
end

for i = 1:2,
    subplot(2, 5, i, 'Spacing', 0.02, 'Margin', 0.02);
    imagesc(D(:,:,i),climst2);
    axis off; % remove the axis ticks and labels
    colormap(T2colormap);
    title('T2')
    daspect([1 1 1])
end

%}

% Set up the 3D array and define the axis and speed
axis = 3;
speed = 5;

% Load the 3D matrix from .mat file
figure(1)
filename = 'AxT1.gif';
check_flag = 0;
for n = 8:22
    subplot();%'Position', [0., 0., 1., 1.]
    data = squeeze(D(160:310,140:290,n));
    imagesc(data,climst1);%axis image;colormap gray;
    axis off;
    colormap(T1colormap);
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
%{
fig = figure();
set(gcf);

v = VideoWriter('CorT2.avi','Uncompressed AVI');
open(v);
NumFrames = 400;

for ii = 100:300,

    subplot();%'Position', [0., 0., 1., 1.]
    data = squeeze(D(ii,:,10:50));
    imagesc(data,climst1);%axis image;colormap gray;
    axis off;
    colormap(T1colormap);
    %title(strcat(num2str(ii),'/',num2str(NumFrames), ' slice'));
    daspect([1 1 1])
    %colorbar
    set(gca, 'box', 'off');    
    drawnow;
    set(gca,'color','white');
    set(gcf,'color','white');
    frame = getframe(fig);
    writeVideo(v,frame);

    
    %im = frame2im(frame);
    %[imind,cm] = rgb2ind(im,256);
    pause(1/speed);
    % Save the frame as a GIF image
    %{
    if ii == 1;
        filename = 'animation.gif';
        imwrite(data,filename,'gif', 'Loopcount',inf, 'DelayTime', 2);
    else
        imwrite(data,filename,'gif','WriteMode','append', 'DelayTime', 2);
    end
    scalebarpsn('location','sw'); hold on;
    %}
end
% Write each frame to the video object

close(v);


% MATLAB program to convert video into slow motion
clc;clear;close all;
  
 % load the video.
obj = VideoReader('AxT2_PF70.avi');  
  
% Write in new variable
obj2= VideoWriter('AxT2_PF70_slow.avi');    
  
% decrease framerate 
obj2.FrameRate = 3;              
open(obj2);
  
% for reading frames one by one
while hasFrame(obj)              
    k = readFrame(obj); 
  
    % write the frames in obj2.         
    obj2.writeVideo(k);          
end
  
close(obj2);
%addpath('./Colormap');





A = t2map3d(100:250,100:250,1);

subplot()
imagesc(A);
title('T2')
colorbar

A2 = 255*(A - min(A(:))) ./ (max(A(:)) - min(A(:))); %scale values between 0 and 255
A2 = cast(A2,'uint8');
imwrite(A2,'T22_deb.png','png'); 
%}
